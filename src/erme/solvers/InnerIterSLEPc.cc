//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   InnerIterSLEPc.cc
 * \author Jeremy Roberts
 * \date   Nov 7, 2011
 * \brief  InnerIterSLEPc member definitions.
 * \note   Copyright (C) 2011 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//
// $Rev::                                               $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date::                                              $:Date of last commit
//---------------------------------------------------------------------------//

#include "serment_config.h"

#ifdef SERMENT_ENABLE_SLEPC

#include <cmath>
#include <iostream>

#include "InnerIterSLEPc.hh"

#include "utilities/DBC.hh"


// Constructor.
InnerIterSLEPc::InnerIterSLEPc(SP_M M, SP_R R)
  : InnerIterBase(M, R)
{
  int ierr;

  // Create the shell matrix.
  Mat MR;
  ierr = MatCreateShell(PETSC_COMM_WORLD, R->m, R->m, R->m, R->m, this, &MR);
  Insist(!ierr, "Error creating MR shell matrix.");
  ierr = MatShellSetOperation(MR, MATOP_MULT, (void(*)(void)) apply_MR);
  Insist(!ierr, "Error setting MR matrix-vector operation.");

  // Create the context.
  ierr = EPSCreate(PETSC_COMM_WORLD, &d_eps);
  Insist(!ierr, "Error creating EPS context.");

  // Set the operator.
  ierr = EPSSetOperators(d_eps, MR, PETSC_NULL);
  Insist(!ierr, "Error setting MR as EPS operator.");

  // Set the problem type.
  ierr = EPSSetProblemType(d_eps, EPS_NHEP);
  Insist(!ierr, "Error setting EPS problem type.");

  // Set the solver type to power and one eigenvalue.
  ierr = EPSSetType(d_eps, EPSKRYLOVSCHUR);
  Insist(!ierr, "Error defaulting EPS to EPSKRYLOVSCHUR.");
  ierr = EPSSetWhichEigenpairs(d_eps, EPS_LARGEST_MAGNITUDE);
  Insist(!ierr, "Error selecting EPS eigenpairs.");

  // Then allow for user choice.
  ierr = EPSSetFromOptions(d_eps);
  Insist(!ierr, "Error setting EPS from options.");

}

// Destructor.
InnerIterSLEPc::~InnerIterSLEPc()
{
  // Nothing here for now.
}

// Solve M*R*J = lambda*J for the dominant eigenpair.
scalar InnerIterSLEPc::solve(int            max_iters,
                             double         tol,
                             SP_vector      J_in,
                             SP_vector      J)
{
  Require(nax_iters > 0);
  Require(tol > 0);
  Require(J_in);            // No null
  Require(J);               // vectors

  scalar lambda_real;
  scalar lambda_imag;
  Vec    J_real;
  Vec    J_imag;
  int    ierr;


  ierr = EPSSetTolerances(d_eps, tol, max_iters);
  Insist(!ierr, "Error setting EPS tolerances.");

  ierr = EPSSetInitialSpace(d_eps, 1, &(J_in->V));
  Insist(!ierr, "Error setting initial guess.");

  ierr = EPSSolve(d_eps);
  Insist(!ierr, "Error solving current eigenvalue problem.");

  int numit = 0;
  ierr = EPSGetIterationNumber(d_eps, &numit);
  Insist(!ierr, "Error getting iteration count.");
  std::cout << "INNERS: Number of iterations =" << numit << std::endl;
  
  // Get the dominant mode.
  MatGetVecs(d_R->M,PETSC_NULL,&J_imag);
  MatGetVecs(d_R->M,PETSC_NULL,&J_real);
  ierr = EPSGetEigenpair(d_eps, 0, &lambda_real, &lambda_imag, J_real, J_imag);
  Insist(!ierr, "Error getting eigenpair.");

  // Copy the dominant mode to the output vector.
  VecCopy(J_real, J->V);
  scalar lambda = lambda_real;

  // Free temporary
  VecDestroy(&J_imag);
  VecDestroy(&J_real);

  return lambda;
}

// Matrix-vector wrapper for M*R
inline PetscErrorCode apply_MR(Mat A, Vec J_in, Vec J_out)
{
  // Get the PETSc context
  void *ctx;
  MatShellGetContext(A, &ctx);
  // and cast it to the base.  This is okay, since we need only apply
  //   M and R, and the base has these.
  InnerIterSLEPc *inner = (InnerIterSLEPc*) ctx;

  // Apply R to X_in and then M to the result.
  MatMult( (inner->R())->M, J_in,                 (inner->J_tmp())->V );
  MatMult( (inner->M())->M, (inner->J_tmp())->V,  J_out               );

  // No errors.
  return 0;
}

#endif

//---------------------------------------------------------------------------//
//              end of InnerIterSLEPc.cc
//---------------------------------------------------------------------------//
