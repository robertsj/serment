//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   EigenSolver.cc
 *  @brief  EigenSolver class definition
 *  @note   Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "EigenSolver.hh"

namespace linear_algebra
{

//----------------------------------------------------------------------------//
EigenSolver::EigenSolver(SP_matrix A,
                         SP_matrix B,
                         int max_iters,
                         double tolerance)
  : d_A(A)
  , d_B(B)
  , d_lambda(0.0)
  , d_lambda_imaginary(0.0)
  , d_number_iterations(0)
  , d_maximum_iterations(max_iters)
  , d_tolerance(tolerance)
{
  // Preconditions
  Require(A);

  // Watch for a generalized problem
  bool generalized = false;
  if (d_B) generalized = true;

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
  ierr = EPSSetWhichEigenpairs(d_solver, EPS_LARGEST_REAL);
  Insist(!ierr, "Error selecting EPS eigenpairs.");

  // Set the tolerances
  ierr = EPSSetTolerances(d_solver, d_tolerance, d_maximum_iterations);
  Insist(!ierr, "Error setting EPS tolerances.");

  // Then allow for user choice.
  ierr = EPSSetFromOptions(d_solver);
  Insist(!ierr, "Error setting EPS from options.");

  // Set temporary storage
  d_x_imaginary = new linear_algebra::Vector(A->number_local_columns());

  // Postconditions
  Ensure(!ierr);
}

//----------------------------------------------------------------------------//
double EigenSolver::solve(SP_vector x)
{
  // Preconditions
  Require(x);

  // Solve the eigenproblem
  PetscErrorCode ierr = EPSSolve(d_solver);
  Insist(!ierr, "Error solving eigenvalue problem.");

  // Extract the number of iterations
  ierr = EPSGetIterationNumber(d_solver, &d_number_iterations);
  Insist(!ierr, "Error getting iteration count.");

  // Check if converged
  int nconv;
  ierr = EPSGetConverged(d_solver, &nconv);
  if (nconv == 0)
  {
    // \todo Need a graceful way to handle this.  It sucks that SLEPc
    //       doesn't allow one to extract the unconverged iterate...
    THROW("SLEPc did not converge!  Exiting...");
  }

  // Get the dominant mode.
  ierr = EPSGetEigenpair(d_solver, 0,
                         &d_lambda, &d_lambda_imaginary,
                         x->V(), d_x_imaginary->V());
  Insist(!ierr, "Error getting eigenpair.");

  // @todo Temporary check for NaN
  volatile double tmp = d_lambda;
  Insist(!(tmp != tmp), "NaN encountered in SLEPc...");
  Insist(d_lambda > 0.0, "Negative lambda...");

  // Return the eigenvalue
  return d_lambda;
}

} // end namespace linear_algebra

//----------------------------------------------------------------------------//
//              end of file EigenSolver.cc
//----------------------------------------------------------------------------//

