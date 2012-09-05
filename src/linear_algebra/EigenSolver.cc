//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   EigenSolver.cc
 * \author robertsj
 * \date   Sep 5, 2012
 * \brief  EigenSolver class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#include "EigenSolver.hh"

namespace linear_algebra
{

EigenSolver::EigenSolver(SP_matrix A, SP_matrix B)
  : d_A(A)
  , d_B(B)
{
  // Preconditions
  Require(A);

  bool generalized = false;
  if (d_B) generalize = true;

  // Create the context.
  PetscErrorCode ierr = EPSCreate(PETSC_COMM_WORLD, &d_solver);
  Insist(!ierr, "Error creating EPS context.");

  // Set the operator.
  if (generalized)
    ierr = EPSSetOperators(d_solver, d_A->A(), d_B->A());
  else
    ierr = EPSSetOperators(d_solver, d_A->A(), PETSC_NULL);
  Insist(!ierr, "Error setting EPS operator.");

  // Set the problem type.
  ierr = EPSSetProblemType(d_solver, EPS_NHEP);
  Insist(!ierr, "Error setting EPS problem type.");

  // Set the solver type to krylovschur and one largest eigenvalue.
  ierr = EPSSetType(d_solver, EPSKRYLOVSCHUR);
  Insist(!ierr, "Error defaulting EPS to EPSKRYLOVSCHUR.");
  ierr = EPSSetWhichEigenpairs(d_solver, EPS_LARGEST_MAGNITUDE);
  Insist(!ierr, "Error selecting EPS eigenpairs.");

  // Then allow for user choice.
  ierr = EPSSetFromOptions(d_solver);
  Insist(!ierr, "Error setting EPS from options.");

  // Postconditions
  Ensure(!ierr);
}

void EigenSolver::solve(SP_vector b, SP_vector x)
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
//              end of file EigenSolver.cc
//---------------------------------------------------------------------------//

