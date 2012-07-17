//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   TwoGroupKernelSolver.cc
 * \author Jeremy Roberts
 * \date   10/14/2011
 * \brief  TwoGroupKernel class member definitions.
 * \note   Copyright (C) 2011 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 140                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-09-14 12:53:40 -0400 (Wed, 14 Sep 2011) $:Date of last commit
//---------------------------------------------------------------------------//



#include "TwoGroupKernelSolver.hh"
#include "Definitions.hh"

#include "utilities/DBC.hh"
#include "utilities/Constants.hh"


namespace erme
{

TwoGroupKernelSolver::TwoGroupKernelSolver(int nh, int order)
  : d_nh(nh), d_n(nh * nh * 2), d_rowindex(d_n + 1, 0), d_rowcount(d_n, 0)
{
  // Initialize rows and columns

  //---------------------------------------------------------------------------//
  //  Generate the row pointers.
  //---------------------------------------------------------------------------//

  getrowcount();  // count = # entries in a row
  getrowindex();  // index tells how many entries preceed this row

  //---------------------------------------------------------------------------//
  //  Generate the column indices.
  //---------------------------------------------------------------------------//
  int nnz = d_rowindex[d_n];               // number of nonzeros
  d_column.resize(nnz);                    // allocate the column indices
  columns(d_column, d_rowindex);           // and compute them.  these are reused.

  //---------------------------------------------------------------------------//
  //  Generate the vectors.
  //---------------------------------------------------------------------------//
  VecCreate( PETSC_COMM_SELF, &v_f);
  VecSetSizes( v_f, PETSC_DECIDE, d_n );
  VecSetFromOptions( v_f );
  VecSet( v_f, 0.0 );
  VecAssemblyBegin( v_f );
  VecAssemblyEnd( v_f );
  VecDuplicate( v_f, &v_si );
  VecDuplicate( v_f, &v_sf );
  VecDuplicate( v_f, &v_b );
  VecDuplicate( v_f, &v_phi);
  VecDuplicate( v_f, &v_phi_old );
  VecDuplicate( v_f, &v_err );

  // Build legendre polynomial for spatial expansion.
  d_legendre.buildMe(d_nh, d_nh, order);

  //---------------------------------------------------------------------------//
  //  Generate and set the matrix.
  //---------------------------------------------------------------------------//
  MatCreateSeqAIJ( PETSC_COMM_SELF, d_n, d_n, 0, &d_rowcount[0], &A );
  MatSeqAIJSetColumnIndices(A, &d_column[0]);
  KSPCreate( PETSC_COMM_SELF, &ksp );
  KSPSetTolerances(ksp,1e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT );
  KSPSetType(ksp,KSPGMRES);
  KSPGMRESSetRestart(ksp,30);
  KSPSetFromOptions(ksp);
  for (int row = 0; row < d_n; row++)
  {
    int rownum = 1;
    int rowidx[] = { row };
    int colnum = d_rowcount[row];
    int colidx[colnum];
    double val[colnum];
    for (int ci = 0; ci < colnum; ci++)
    {
      colidx[ci] = d_column[d_rowindex[row] + ci];
      val[ci]    = 1.0; // set to ones to ensure right nonzero structure.
    }
    MatSetValues(A, rownum, rowidx, colnum, colidx, val, INSERT_VALUES);
  }
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  KSPSetOperators( ksp, A, A, SAME_NONZERO_PATTERN );
}

void TwoGroupKernelSolver::update_matrix(Vec_Dbl &data)
{
  MatGetArray(A, &a_A);
  // Reset to zero.
  for (int i = 0; i < d_nnz; i++)
    a_A[i] = 0.0;
  VecGetArray(v_f, &a_f);
  matrix(data, a_A, a_f);
  VecRestoreArray(v_f, &a_f);
  MatRestoreArray(A, &a_A);
}

void TwoGroupKernelSolver::set_source(Vec_Dbl &data, int group, int order)
{
  // We set only the data on the bottom edge.
  int nh = d_nh;
  double dim    = data[NODEDIM];  // node length
  double D[]    =
  { data[D1], data[D2] };         // group diffusion coefficients
  double h = dim / d_nh;          // mesh size

  // Get arrays
  double *a_legendre;
  VecGetArray(d_legendre.Px[order], &a_legendre);
  VecGetArray(v_si, &a_si);

  // Set the boundary for this group and wipe the source for the other group
  int row = 0;
  for (int g = 0; g < 2; g++)
  {
    for (int i = 0; i < d_nh; i++)
    {
      row = i + g*nh*nh;
      if (g == group)
        a_si[row] = 8.0 * D[g] * h * a_legendre[i] / (4.0 * D[g] + h);
      else
        a_si[row] = 0.0;
    }
  }
  // Replace arrays
  VecRestoreArray(d_legendre.Py[order], &a_legendre);
  VecRestoreArray(v_si, &a_si);

}

void TwoGroupKernelSolver::solve(int o, int g)
{

  KSPSetOperators(ksp, A, A, SAME_NONZERO_PATTERN);

  double volume = inputs[state][1];
  volume *= volume;

  // Reset fission source and right hand side.  Note, 1/k is included in operator.
  VecCopy(v_f, v_sf); // sf is cross section; should help a bit;
  VecSet(v_b, 0.0);   // zeros

  VecAssemblyBegin(v_si);
  VecAssemblyEnd(v_si);
  VecAssemblyBegin(v_sf);
  VecAssemblyEnd(v_sf);
  VecAssemblyBegin(v_b);
  VecAssemblyEnd(v_b);

  int    i_norm_err = 1;
  double norm_err   = 1.0;
  bool   converged  = false;

  // Now iterate on fission source
  for (int it = 0; it < 100; it++)
  {

    // Build right hand side
    VecCopy(v_sf, v_b);      // b <- sf
    VecAXPY(v_b, 1.0, v_si); // b <- b + si

    // Solve the problem
    KSPSolve(ksp, v_b, v_phi);

    // compute new fission source
    VecGetArray(v_sf,  &a_sf);
    VecGetArray(v_phi, &a_phi);
    VecGetArray(v_f,   &a_f);
    for (int i = 0; i < n / 2; i++)
      a_sf[i] = a_f[i] * a_phi[i]; // group 1
    for (int i = 0; i < n / 2; i++)
      a_sf[i] += a_f[i + n / 2] * a_phi[i + n / 2]; // group 2
    VecRestoreArray(v_sf,  &a_sf);
    VecRestoreArray(v_phi, &a_phi);
    VecRestoreArray(v_f,   &a_f);

    // check phi convergence
    VecCopy(v_phi_old, v_err);         // err <- phi_old
    VecAXPY(v_err, -1.0, v_phi);       // err <- phi_old - phi
    VecNorm(v_err, NORM_1, &norm_err); // ||err||_oo (= max(|err[i]|,i=1,n)
    VecCopy(v_phi, v_phi_old);
    if (norm_err < 1.0e-4)
    {
      //VecView(v_phi, PETSC_VIEWER_STDOUT_SELF);
      //cout << "converged at it " << it << endl;
      converged = true;
      break;
    }
  } // end fission loop
  if (!converged)
  {
    cout << " warning: fission didn't converge for state " << state << endl;
  }

  //-----------------------------------------------------------------
  // FISSION RESPONSE = volume-integrated fission source
  double tmp;
  VecSum(v_sf, &tmp);
  output[16 + g] = tmp * volume;

  //-----------------------------------------------------------------
  // ABSORPTION RESPONSE = volume-integrated absorption rate
  VecGetArray(v_phi, &a_phi);
  tmp = 0;
  double siga = inputs[state][4] - inputs[state][8];
  for (int i = 0; i < n / 2; i++)
    tmp += a_phi[i] * siga;
  siga = inputs[state][5];
  for (int i = n / 2; i < n; i++)
    tmp += a_phi[i] * siga;
  output[18 + g] = tmp * volume;

  //-----------------------------------------------------------------
  // FLUX RESPONSE = volume-integrated flux
  double tmp;
  VecSum(v_sf, &tmp);
  output[16 + g] = tmp * volume;


}


// implementation

// data, matrix nonzeros, fission operator,
void TwoGroupKernelSolver::matrix(Vec_Dbl &data, double *t, double *f)
{
  using namespace constants;

  // set data
  int nh = d_nh;
  double keff   = data[KEFF];       // user set keff
  double dim    = data[NODEDIM];    // node length
  double D[]    =
  { data[D1], data[D2] };           // group diffusion coeficients
  double R[]    =
  { data[R1], data[A2] };           // group removal cross sections
  double F[]    =
  { data[F1]/keff, data[F2]/keff }; // group nu * fission cross sections / keff
  double S12    = data[8];      // g1 to g2 scattering cross section
  int nh_sq     = nh * nh;      //
  double h      = dim / nh;     // delta
  int ndiags    = 6;            // number of diagonals (one lower for scatter)
  int n         = 2*nh*nh;      // matrix size (# groups x # cells)
  double ntwo   = -2.0;         //
  double h_sq   = h * h;

  int k; // this is the row index
  int j = 0;
  int i = 0;
  int g = 0;
  int kgp = 0;
  double n_two_D = 0;

  for (int row = 0; row < n; row++)
  {
    g       = row/(nh*nh); // = 0 if row < nh_sq, 1 otherwise
    i       = row % nh;
    j       = ((row - i)/nh) % nh;
    n_two_D = ntwo * D[g]; // -2*diffusion
    k       = d_rowindex[row];
    kgp     = k + g;

    if (g == 1)
    {
      t[k] = -h_sq * S12;
    }
    f[row]   = h_sq * F[g] / keff; // fission goes everywhere
    // interior coefficients [B-0,L-1,C-2,R-3,T-4]
    if (i > 0 and i < nh - 1 and j > 0 and j < nh - 1)
    {
      t[kgp + 0] = t[kgp + 1] = t[kgp + 3] = t[kgp + 4] = -D[g];
      t[kgp + 2] = h_sq * R[g] + 4 * D[g];
    }
    // left edge             [B-0,C-1,R-2,T-3]
    else if (i == 0 and j > 0 and j < nh - 1)
    {
      t[kgp + 0] = t[kgp + 2] = t[kgp + 3] =  -D[g];
      t[kgp + 1] = h_sq * R[g] + D[g] * (3 + h * 2 / (4 * D[g] + h));
    }
    // right edge            [B-0,L-1,C-2,T-3]
    else if (i == nh - 1 and j > 0 and j < nh - 1)
    {
      t[kgp + 0] = t[kgp + 1] = t[kgp + 3] =  -D[g];
      t[kgp + 2] = h_sq * R[g] + D[g] * (3 + h * 2 / (4 * D[g] + h));
    }
    // bottom edge           [L-0,C-1,R-2,T-3]
    else if (j == 0 and i > 0 and i < nh - 1)
    {
      t[kgp + 0] = t[kgp + 2] = t[kgp + 3] =  -D[g];
      t[kgp + 1] = h_sq * R[g] + D[g] * (3 + h * 2 / (4 * D[g] + h));
    }
    // top edge              [B-0,L-1,C-2,R-3]
    else if (j == nh - 1 and i > 0 and i < nh - 1)
    {
      t[kgp + 0] = t[kgp + 1] = t[kgp + 3] =  -D[g];
      t[kgp + 2] = h_sq * R[g] + D[g] * (3 + h * 2 / (4 * D[g] + h));
    }
    // bottom left corner    [C-0,R-1,T-2]
    else if (i == 0 and j == 0)
    {
      t[kgp + 1] = t[kgp + 2] =  -D[g];
      t[kgp + 0] = h_sq * R[g] + 2 * D[g] + D[g] * h * 4 / (4 * D[g] + h);
    }
    // top left corner       [B-0,C-1,R-2]
    else if (i == 0 and j == nh - 1)
    {
      t[kgp + 0] = t[kgp + 2] =  -D[g];
      t[kgp + 1] = h_sq * R[g] + 2 * D[g] + D[g] * h * 4 / (4 * D[g] + h);
    }
    // bottom right corner   [L-0,C-1,T-2]
    else if (i == nh - 1 and j == 0)
    {
      t[kgp + 0] = t[kgp + 2] =  -D[g];
      t[kgp + 1] = h_sq * R[g] + 2 * D[g] + D[g] * h * 4 / (4 * D[g] + h);
    }
    // top right corner      [B-0,L-1,C-2]
    else if (i == nh - 1 and j == nh - 1)
    {
      t[kgp + 0] = t[kgp + 1] =  -D[g];
      t[kgp + 2] = h_sq * R[g] + 2 * D[g] + D[g] * h * 4 / (4 * D[g] + h);
    }
    else
    {

    }
  }
}

void TwoGroupKernelSolver::columns(Vec_Int c, Vec_Int rp)
{
  int nh = d_nh;
  int nh_sq = nh * nh;
  int n     = 2*nh_sq;
  int k     = 0;
  int j     = 0;
  int i     = 0;
  int g     = 0;
  int kgp   = 0;

  for (int row = 0; row < n; row++)
  {
    if (row >= n/2)
      g = 1;
    i = row % nh;
    j = ((row - i) / nh) % nh;
    k = rp[row];
    kgp = k + g;
    if (g == 1)
      c[k] = row - nh * nh;

    // interior coefficients [B-0,L-1,C-2,R-3,T-4]
    if (i > 0 and i < nh - 1 and j > 0 and j < nh - 1)
    {
      c[kgp + 0] = row - nh;
      c[kgp + 1] = row - 1;
      c[kgp + 2] = row;
      c[kgp + 3] = row + 1;
      c[kgp + 4] = row + nh;
    }
    // left edge
    else if (i == 0 and j > 0 and j < nh - 1)
    {
      c[kgp + 0] = row - nh;
      c[kgp + 1] = row;
      c[kgp + 2] = row + 1;
      c[kgp + 3] = row + nh;
    }
    // right edge
    else if (i == nh - 1 and j > 0 and j < nh - 1)
    {
      c[kgp + 0] = row - nh;
      c[kgp + 1] = row - 1;
      c[kgp + 2] = row;
      c[kgp + 3] = row + nh;
    }
    // bottom edge
    else if (j == 0 and i > 0 and i < nh - 1)
    {
      c[kgp + 0] = row - 1;
      c[kgp + 1] = row;
      c[kgp + 2] = row + 1;
      c[kgp + 3] = row + nh;
    }
    // top edge
    else if (j == nh - 1 and i > 0 and i < nh - 1)
    {
      c[kgp + 0] = row - nh;
      c[kgp + 1] = row - 1;
      c[kgp + 2] = row;
      c[kgp + 3] = row + 1;
    }
    // bottom left corner
    else if (i == 0 and j == 0)
    {
      c[kgp + 0] = row;
      c[kgp + 1] = row + 1;
      c[kgp + 2] = row + nh;
    }
    // top left corner
    else if (i == 0 and j == nh - 1)
    {
      c[kgp + 0] = row - nh;
      c[kgp + 1] = row;
      c[kgp + 2] = row + 1;
    }
    // bottom right corner [L-0,C-1,T-2]
    else if (i == nh - 1 and j == 0)
    {
      c[kgp + 0] = row - 1;
      c[kgp + 1] = row;
      c[kgp + 2] = row + nh;
    }
    // top right corner [B-0,L-1,C-2]
    else if (i == nh - 1 and j == nh - 1)
    {
      c[kgp + 0] = row - nh;
      c[kgp + 1] = row - 1;
      c[kgp + 2] = row;
    }
    else
    {

    }
  }
}

int TwoGroupKernelSolver::countrow(int row)
{
  // indices
  int g = 0; // = 0 if row < nh_sq, 1 otherwise
  if (row >= d_nh*d_nh)
    g = 1;
  int i = row % d_nh;
  int j = ((row - i)/d_nh) % d_nh;

  // count number in my row
  int count = -10;
  if (i > 0 and i < d_nh - 1 and j > 0 and j < d_nh - 1)            // internal
    count = 5;
  else if ((i == 0 or i == d_nh - 1) and j > 0 and j < d_nh - 1)    // edge
    count = 4;
  else if ((j == 0 or j == d_nh - 1) and i > 0 and i < d_nh - 1)    // edge
    count = 4;
  else if ((i == 0 or i == d_nh - 1) and (j == 0 or j == d_nh - 1)) // corner
    count = 3;
  if (g == 1) count++;
  return count;
}

void TwoGroupKernelSolver::getrowcount()
{
  // These will remain constant for all problems.
  for (int i = 0; i < d_n; i++)
    d_rowcount[i] = countrow(i);
}

void TwoGroupKernelSolver::getrowindex()
{
  // These will remain constant for all problems.
  d_rowindex[0] = 0;
  for (int i = 0; i < d_n; i++)
    d_rowindex[i+1] = d_rowindex[i] + countrow(i);
}


} // end namespace erme

//---------------------------------------------------------------------------//
//                 end of TwoGroupKernel.hh
//---------------------------------------------------------------------------//
