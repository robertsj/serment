//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LinearSolver.cc
 * \author robertsj
 * \date   Sep 5, 2012
 * \brief  LinearSolver class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#include "LinearSolver.hh"

namespace linear_algebra
{

LinearSolver::LinearSolver(SP_matrix A, SP_matrix P)
  : d_A(A)
  , d_P(P)
{
  // Preconditions
  Require(A);
  Require(P);

  // Create the KSP object.
  PetscErrorCode ierr = KSPCreate(PETSC_COMM_WORLD, &d_solver);
  Insist(!ierr, "Error creating KSP solver.");

  // Set the operator.
  ierr = KSPSetOperators(d_solver, A->A(), P->A(), SAME_NONZERO_PATTERN);
  Insist(!ierr, "Error setting KSP operators.");

  // Set tolerances.
  ierr = KSPSetTolerances(d_solver, PETSC_DEFAULT,
                          PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);

  // Allow for command line flags.
  ierr = KSPSetFromOptions(d_solver);

  // Postconditions
  Ensure(!ierr);
}

void LinearSolver::solve(SP_vector b, SP_vector x)
{
  // Preconditions
  Require(b);
  Require(x);

  // Solve
  PetscErrorCode ierr = KSPSolve(d_solver, b->V(), x->V());

  // Postconditions
  Insist(!ierr, "Error in KSPSolve.");
}

} // end namespace linear_algebra

//---------------------------------------------------------------------------//
//              end of file LinearSolver.cc
//---------------------------------------------------------------------------//

