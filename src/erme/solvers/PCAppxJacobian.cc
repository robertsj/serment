//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PCAppxJacobian.cc
 * \author Jeremy Roberts
 * \date   10/25/2010
 * \brief  Member definitions of abstract class PCAppxJacobian
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 177                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-12-09 16:31:49 -0500 (Fri, 09 Dec 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#include <iostream>
#include "LinAlg.hh"
#include "GlobalProblem.hh"
#include "Newton.hh"
#include "PCAppxJacobian.hh"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//


PCAppxJacobian::PCAppxJacobian( integer        a,
                                integer        b,
                                integer        c,
                                integer        const nnz[],
                                SermentVector *unk,
                                GlobalProblem *pr )
  : SermentMatrixCRS(a,b,c,nnz),
    x(unk),
    problem(pr)
{
    // The Jacobian has the form
  //
  //   | (M*R-lambda*I)   M*R_k*Jinc           -Jinc |
  //   | (F - k*L)        (F_k-k*L_k-L)*Jinc    0    |
  //   | -Jinc'           0                     0    |
  //

  // get the unknowns into an array
  scalar *x_a;
  VecGetArray( x->V, &x_a );
  scalar k = x_a[m - 2]; // keff
  scalar lambda = x_a[m - 1]; // current eigenvalue

  problem->R->updateData(k);
  problem->L.updateData(k);
  problem->F.updateData(k);
  problem->A.updateData(k);

  scalar time = MPI_Wtime();
  // To build the approximate Jacobian, we need to compute MR = M*R
  //   Since R is stored as BCSR, we need to convert it into CSR.  I *assume*
  //   this cost is small relative to the benefit of using R in block
  //   format.
  Mat MR;
  Mat Rcsr;
  MatConvert(problem->R->M, MATAIJ, MAT_INITIAL_MATRIX, &Rcsr);

  time = MPI_Wtime() - time;
  cout << " time convert = " << time << endl;

  time = MPI_Wtime();

  PetscViewer bview;
  if (true)
  {
    MatMatMult(problem->M->M, Rcsr, MAT_INITIAL_MATRIX, 1.0, &MR);
    //        PetscViewerBinaryOpen(PETSC_COMM_WORLD,"MR->bin",FILE_MODE_WRITE,&bview);
    //        MatView(MR,bview);
  }
  else
  {
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "MR->bin", FILE_MODE_READ, &bview);
    //MatLoad(bview, MATSEQAIJ, &MR);
  }
  time = MPI_Wtime() - time;
  cout << " time MR = " << time << endl;

  time = MPI_Wtime();
  // now create M*R-lambda*I
  SermentVector D1(m - 2);
  for (integer i = 0; i < m - 2; i++)
    D1.insertVal(i, -lambda);
  // MatDiagonalSet( MR, D1.V, ADD_VALUES );

  SermentVector D2(m - 2);
  MatGetDiagonal(MR, D2.V);
  //    cout << " D2 " << endl;
  //    D2.viewMe();
  D2.vecAYPV(1.0, D1);
  scalar *d_a;
  VecGetArray( D2.V, &d_a );
  //    cout << " D2 " << endl;
  //    D2.viewMe();


  //  MatView( MR, PETSC_VIEWER_STDOUT_SELF );

  // destroy Rcsr, since we don't need it anymore
  MatDestroy(&Rcsr);
  // MatDestroy( Mbcsr );
  time = MPI_Wtime() - time;
  cout << " diagonal = " << time << endl;

  time = MPI_Wtime();

  // now we begin building the main block of the PC matrix going row by row
  // of MR.
  integer row, ncols;
  const integer *cols;
  const scalar *vals;
  scalar tmpval = 0.0;
  for (row = 0; row < m - 2; row++)
  {
    MatGetRow(MR, row, &ncols, &cols, &vals);
    //this->insertVals(  vals, 1, row, ncols, cols );
    //        for ( int j = 0; j < ncols; j++ )
    //            cout << " row = " << row << " col = " << cols[j] << " v = " << vals[j] << endl;
    MatSetValues(M, 1, &row, ncols, cols, vals, INSERT_VALUES);
    tmpval = d_a[row];
    MatSetValues(M, 1, &row, 1, &row, &tmpval, INSERT_VALUES); // just insert it over a possibly inserted values of MR[i,i]
    MatRestoreRow(MR, row, &ncols, &cols, &vals);
  }
  MatDestroy(&MR);
  VecRestoreArray( D2.V, &d_a );
  D1.releaseMe();
  D2.releaseMe();

  time = MPI_Wtime() - time;
  cout << " time get rows = " << time << endl;
  time = MPI_Wtime();

  // now insert the -J terms
  integer *idx;
  idx = new integer[m - 2];
  idx[0] = 0;
  for (integer i = 1; i < m - 2; i++)
    idx[i] = idx[i - 1] + 1;
  SermentVector tmp1(m - 2);
  SermentVector tmp2(m - 2);

  VecPlaceArray(tmp1.V, x_a); // J
  //    cout << "J=" << endl;
  //    tmp1.viewMe();
  VecSet(tmp2.V, 0.0);
  tmp2.vecAYPV(-1.0, tmp1); // tmp2 = -J

  scalar *tmp_a;
  VecGetArray( tmp2.V, &tmp_a ); // put -J in matrix
  row = m - 1;
  this->insertVals(tmp_a, 1, &row, m - 2, idx);
  this->insertVals(tmp_a, m - 2, idx, 1, &row);
  VecRestoreArray( tmp2.V, &tmp_a );

  // compute and insert the middle term, -(Abs+Leak)J
  scalar absleak =
      -(problem->A.vecDot(tmp1) + problem->L.computeLeakage(tmp1));
  this->insertVal(absleak, m - 2, m - 2);
  VecResetArray(tmp1.V);

  // finally, compute (F - k*(abs+leak) );
  // get the leakage vector (not stored explicitly)
  problem->L.getLeakageVec(tmp1);  // = L
  tmp1.vecAYPV(1.0, problem->A);   // = L + A
  tmp1.vecScale(-k);               // = -k(L+A)
  tmp1.vecAYPV(1.0, problem->F);   // = F - k(L+A)
  VecGetArray( tmp1.V, &tmp_a );
  row = m - 2;
  this->insertVals(tmp_a, 1, &row, m - 2, idx);
  VecRestoreArray( tmp1.V, &tmp_a );

  // put a one in the bottom right hand corner or petsc complains
  row = m - 1;
  this->insertVal(0.00000001, m - 1, m - 1);

  // get rid of the two temporary vectors
  tmp1.releaseMe();
  tmp2.releaseMe();

  time = MPI_Wtime() - time;
  cout << " time for rest = " << time << endl;
  time = MPI_Wtime();

  // we forgo (for now) any finite differencing and set the matrix
  this->checkReady();

    // debug
//    cout << " i'm approx: " << endl;
//    viewMe();
//    PetscViewer viewer;
//    PetscViewerASCIIOpen(PETSC_COMM_WORLD,"jacappx.output",&viewer);
//    PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_DENSE );
//    MatView(M,viewer);
}

//---------------------------------------------------------------------------//
//                 end of PCAppxJacobian.cc
//---------------------------------------------------------------------------//

