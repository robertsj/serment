//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Newton.cc
 * \author Jeremy Roberts
 * \date   11/24/2010
 * \brief  Member definitions of base class Newton
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 178                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-12-12 13:36:39 -0500 (Mon, 12 Dec 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#include <iostream>
#include <cmath>
#include "petscsnes.h"
#include "LinAlg.hh"
#include "GlobalInput.hh"
#include "ResponseFunctionServer.hh"
#include "ResponseMatrix.hh"
#include "ResponseMatrixFull.hh"
#include "AbsorptionResponse.hh"
#include "FissionResponse.hh"
#include "LeakageResponse.hh"
#include "ConnectMatrix.hh"
#include "Connect2dCart.hh"
#include "GlobalProblem.hh"
#include "GlobalSolver.hh"
#include "Newton.hh"
#include "ResidualWrap.hh"
#include "JacobianShell.hh"
#include "JacobianEmpty.hh"
#include "PCAppxJacobian.hh"
//#include "PIPCShell.hh"
#include "utilities/DBC.hh"


//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//


Newton::Newton(SP_globalproblem problem, SP_globalinput input) :
  GlobalSolver(problem, input)
{
  std::cout << " CONSTRUCTING Newton " << std::endl;
  // nothing more here for now
}

//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
Newton::~Newton()
{
  // nothing more here right now, but should write Newton->destroy();
  return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief This is the Newton solver.
 *
 */
void Newton::solve()
{

  scalar one = 1.0;
  scalar zero = 0.0;

  // - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
  //  Create nonlinear solver context
  // - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
  SNESCreate(PETSC_COMM_WORLD, &snes);

  // - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
  //  Create vectors for solution and nonlinear function
  // - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
  x = new SermentVector(d_input->degfree + 2);
  f = new SermentVector(d_input->degfree + 2);
  // initialize
  for (integer i = 0; i < d_input->degfree + 2; i++)
  {
    x->insertVal(i, zero);
    f->insertVal(i, one);
  }
  // uniform zeroth order initial guess
  for (integer i = 0; i < d_input->degfree; i = i + (d_input->spaceord + 1)
      * (d_input->angleord + 1) * d_input->numgroups)
  {
    x->insertVal(i, one);
  }
  x->vecScale(1.0 / sqrt(x->vecDot(*x))); // norm'd uniform 0-order guess
  x->insertVal(d_input->degfree, d_input->keff); // keff guess from input
  x->insertVal(d_input->degfree + 1, one); // guess for lambda is unity

  // perform one crude power iteration as initial seed
  scalar keff = d_input->keff;
  scalar lambda = 0;
  scalar *x_a;
  VecGetArray( x->V, &x_a );
  SermentVector Jinc(d_input->degfree);
  SermentVector temp(d_input->degfree);
  VecPlaceArray(Jinc.V, x_a);
  for (int powit = 0; powit < 1; powit++)
  {
    d_problem->R->updateData(keff);
    d_problem->L.updateData(keff);
    d_problem->A.updateData(keff);
    d_problem->F.updateData(keff);
    for (integer i = 0; i < 1000; i++)
    {
      d_problem->R->matVec(Jinc, temp);
      d_problem->M->matVec(temp, Jinc);
      lambda = sqrt(Jinc.vecDot(Jinc));
      Jinc.vecScale(1.0 / lambda);
    }
    keff = d_problem->F.vecDot(Jinc) / (d_problem->A.vecDot(Jinc)
        + d_problem->L.computeLeakage(Jinc));
  }
  VecResetArray(Jinc.V);
  x_a[d_input->degfree] = keff;
  x_a[d_input->degfree + 1] = lambda;
  VecRestoreArray( x->V, &x_a );
  Jinc.releaseMe();
  temp.releaseMe();

  // - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
  //  Set the residual and vector in SNES
  // - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
  SNESSetFunction(snes, f->V, ResidualWrap, this);

  // - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
  //  Set the solver context and solve the problem
  // - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
  KSP ksp;
  PC pc;
  SNESGetKSP(snes, &ksp);
  KSPSetType(ksp, KSPGMRES);

  // - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
  //  Set the preconditioner
  // - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - -

  KSPGetPC(ksp, &pc);

  PCAppxJacobian *JacPC;
  //    PIPCShell *PIPC;


  if (d_input->pctype == 0) // no preconditioner
  {
    PCSetType(pc, PCNONE);
  }
  else if (d_input->pctype == 1) // build the approximate jacobian matrix
  {
    // the approximate jacobian is sparse except for the final two rows
    // and columns.  hence, it makes sense to allocate using the nzz
    // rather than nz parameter.  For the main block, each row has at most
    // degfree/numel + 2 nonzeros;  For the m-1th, m-1, and mth, m-2.
    integer nzz[d_input->degfree + 2];
    for (integer i = 0; i < d_input->degfree; i++)
      nzz[i] = d_input->degfree / d_input->numel + 2;
    nzz[d_input->degfree] = d_input->degfree + 1;
    nzz[d_input->degfree + 1] = d_input->degfree;

    //        for (integer i=0; i < d_input->degfree+2; i++)
    //            cout << " NZZ = " << nzz[i] << endl;
    const integer *nzzp = nzz;
    JacPC = new PCAppxJacobian(d_input->degfree + 2, d_input->degfree + 2, 3,
                               nzzp, x, &(*d_problem));
    SNESSetLagPreconditioner(snes, -1); // -1 indicates we never rebuild
    PCSetType(pc, PCILU);
    PCFactorSetLevels(pc, d_input->ilulevel);
  }
  else if (d_input->pctype == 2) // PI preconditioner
  {

    //        cout << " pc type 2 " << endl;
    //        PIPC = new PIPCShell(d_input->degfree+2,
    //                             d_input->degfree+2,PETSC_NULL,x,f,problem,snes);
    //        PCSetType(pc,PCSHELL);
    //        PCShellSetContext(pc,PIPC);
    //        PCShellSetApply(pc,PIPCMatVecWrap);
    //   //     KSPSetPreconditionerSide(ksp,PC_RIGHT);
  }

  // - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
  //  Create Jacobian matrix data structure
  // - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
  Mat Jac;
  integer jj = 0;
  if (jj == 0)
  {
    fp = new JacobianShell(d_input->degfree + 2, d_input->degfree + 2, this,
                           x, f, &(*d_problem), snes);
    if (d_input->pctype == 0) SNESSetJacobian(snes, fp->M, fp->M,
                                              JacobianEmpty, this);
    else if (d_input->pctype == 1) SNESSetJacobian(snes, fp->M, JacPC->M,
                                                   JacobianEmpty, this);
    else if (d_input->pctype == 2) Insist(d_input->pctype < 2,
                                          "PCTYPE 2 is not implemented yet.");
    //SNESSetJacobian( snes, fp->M, PIPC->M, JacobianEmpty, this );
  }
  else if (jj == 1)
  {
    // testing with matrix free finite difference
    MatCreateSNESMF(snes, &Jac);
    MatMFFDSetFunction(
                       Jac,
                       (PetscErrorCode(*)(void*, Vec, Vec)) SNESComputeFunction,
                       snes);
    SNESSetJacobian(snes, Jac, JacPC->M, MatMFFDComputeJacobian, this);
  }
  else
  {
    // testing with full finite difference jacobian
    MatCreate(PETSC_COMM_SELF, &Jac);
    MatSetSizes(Jac, PETSC_DECIDE, PETSC_DECIDE, d_input->degfree + 2,
                d_input->degfree + 2);
    MatSetFromOptions(Jac);
    SNESSetJacobian(snes, Jac, Jac, SNESDefaultComputeJacobian, this);
  }

  // use the following to quit after just one A*v.  Usually quits newton, too
  //KSPSetTolerances( ksp, 1.e-5, PETSC_DEFAULT, PETSC_DEFAULT, 1 );

  integer its, itslin;
  scalar tol = 1e-8;
  SNESSetTolerances(snes, tol, tol, tol, 2, 1000);
  SNESSetFromOptions(snes);
  SNESLineSearchSet(snes, SNESLineSearchNo, PETSC_NULL);
  SNESSolve(snes, PETSC_NULL, x->V);

  SNESGetIterationNumber(snes, &its);
  SNESGetLinearSolveIterations(snes, &itslin);

  // test fissionrates
//  SermentVector Jo(x->Length() - 2);
//  scalar *j_a;
  VecGetArray( x->V, &x_a ); // get the array
//  VecGetArray( Jo.V, &j_a ); // get the array
//  for (int i = 0; i < Jo.Length(); i++)
//    j_a[i] = x_a[i];
//  VecRestoreArray( x->V, &x_a ); // get the array
//  VecRestoreArray( Jo.V, &j_a ); // get the array
 // fissionRates(&Jo);

  cout << " FINAL NEWTON EIGENVALUES: " << endl;
  printf(" **** FINAL KEFF   = %12.9f \n", x_a[d_input->degfree]);
  printf(" **** FINAL LAMBDA = %12.9f \n", x_a[d_input->degfree + 1]);
  printf(" **** NONLINEAR ITERATIONS  = %8i \n", its);
  printf(" **** LINEAR ITERATIONS     = %8i \n", itslin);

  if (jj == 2) // debug, prints fd-jacobian to ascii and binary
  {
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "jac.output", &viewer);
    PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_DENSE);
    MatView(Jac, viewer);
    MatView(Jac, PETSC_VIEWER_STDOUT_SELF);
    PetscViewer bview;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "jac.bin", FILE_MODE_WRITE, &bview);
    MatView(Jac, bview);
  }

  return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief This function computes the nonlinear residual.
 *
 * The nonlinear residual is defined as
 *   \f[
 *       \mathbf{f(x)} = \left [\begin{array}{c}
 *	        (\mathbf{M}\mathbf{R}(k)-\lambda \mathbf{I}) \mathbf{J_-} \\
 *	        \mathbf{F}(k)\mathbf{J_-} - (k\mathbf{L}(k)\mathbf{J_-} ) \\
 *	        \frac{1}{2} \mathbf{J^T_-} \mathbf{J_-} - \frac{1}{2}  
 *	      \end{array} 
 *       \right ]  = \mathbf{0} \, ,
 *   \f]
 * which is the same as used in the
 *
 * \todo Residual should be a class so that the temp vectors don't have to
 * be rebuilt every time
 *
 */
PetscErrorCode Newton::Residual(SNES snes, Vec X, Vec F, void *ptr)
{

  scalar *x_a, *f_a;
  VecGetArray( X, &x_a );
  VecGetArray( F, &f_a );

  SermentVector Jinc(d_input->degfree);
  SermentVector F_a(d_input->degfree);
  SermentVector temp(d_input->degfree);

  scalar k, lambda;

  // trying vecplacearray with x_a and f_a into smaller arrays
  VecPlaceArray(Jinc.V, x_a);
  VecPlaceArray(F_a.V, f_a);

  k = x_a[d_input->degfree];
  lambda = x_a[d_input->degfree + 1];

  d_problem->R->updateData(k);
  d_problem->L.updateData(k);
  d_problem->F.updateData(k);
  d_problem->A.updateData(k);

  // first m-2 entries
  d_problem->R->matVec(Jinc, temp);
  d_problem->M->matVec(temp, F_a);
  F_a.vecAYPV(-lambda, Jinc);

  // (m-1)th entry
  f_a[d_input->degfree] = d_problem->F.vecDot(Jinc) - k * (d_problem->A.vecDot(Jinc)
      + d_problem->L.computeLeakage(Jinc));
  // (m)th entry
  f_a[d_input->degfree + 1] = 0.5 - 0.5 * Jinc.vecDot(Jinc);

  VecResetArray(Jinc.V);
  VecResetArray(F_a.V);

  VecRestoreArray( X, &x_a );
  VecRestoreArray( F, &f_a );

  Jinc.releaseMe();
  F_a.releaseMe();
  temp.releaseMe();

  return 0;
}

//---------------------------------------------------------------------------//
//                 end of Newton.cc
//---------------------------------------------------------------------------//

