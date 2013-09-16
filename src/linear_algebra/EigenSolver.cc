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
                         double tolerance,
                         const std::string type)
  : d_type(type)
  , d_A(A)
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

  if (d_type != "simplepower")
  {

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
    ierr = EPSSetType(d_solver, type.c_str());
    Insist(!ierr, "Error setting EPS to"+type);
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
  else
  {
    Insist(!d_B, "Simple power is not for generalized problems.");

    // Use this for temporary
    d_x_imaginary = new linear_algebra::Vector(A->number_local_columns());
  }

}

//----------------------------------------------------------------------------//
EigenSolver::~EigenSolver()
{
  if (d_type != "simplepower")
  {
    EPSDestroy(&d_solver);
  }
}

//----------------------------------------------------------------------------//
double EigenSolver::solve(SP_vector x)
{
  // Preconditions
  Require(x);

  if (d_type != "simplepower")
  {
    // Set the initial guess
    Vec V[] = {x->V()};
    PetscErrorCode ierr = EPSSetInitialSpace(d_solver, 1, V);

    // Solve the eigenproblem
    ierr = EPSSolve(d_solver);
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
    //


    Insist(!ierr, "Error getting eigenpair.");

    // @todo Temporary check for NaN
    volatile double tmp = d_lambda;
    Insist(!(tmp != tmp), "NaN encountered in SLEPc...");
    Insist(d_lambda > 0.0, "Negative lambda...");

  }
  else
  {

    PetscBool do_monitor = PETSC_FALSE, has_shift = PETSC_FALSE;
    int ierr = PetscOptionsHasName("", "-eps_monitor", &do_monitor);
    double shift = 1.0;
    ierr = PetscOptionsGetReal("", "-power_shift", &shift, &has_shift);

    SP_vector x0 = d_x_imaginary;
    SP_vector swap;

    // Eigenvector errors.  We'll estimate the dominance ratio using
    // a few iterations, and then estimate how many iterations until
    // the error has died out such that the tolerance is met.
    double x_err = 0.0;

    // Normalize guess to unity
    d_lambda = x->norm(x->L2);
    x->scale(1.0 / d_lambda);

    x0->copy(*x);

    // Do power iterations
    d_number_iterations = 0;
    for (; d_number_iterations < d_maximum_iterations; ++d_number_iterations)
    {
      // compute x <-- (I + A) * x0
      d_A->multiply(*x0, *x);
      if (shift != 0.0)
        x->add_a_times_x(shift, *x0);

      // Compute the eigenvalue
      d_lambda = x->norm(x->L2);

      // Compute residual norm
      x_err = x0->norm_residual(*x, x->L2);

      // Normalize and swap
      x->scale(1.0 / d_lambda);
      x_err = x0->norm_residual(*x, x->L2);
      x0->copy(*x);
      //x->swap(x0);

      // Check for convergence.
      if (do_monitor == PETSC_TRUE)
      {
        printf("simplepower   %6i %20.12e %20.12e \n",
               d_number_iterations, d_lambda-shift, x_err);
      }
      if (x_err < d_tolerance) break;
    }

    d_lambda -= shift;
//    if (d_number_iterations % 2)
//      swap_vector(x.bp(), x0.bp(), swap.bp());

  }

  // Return the eigenvalue
  return d_lambda;
}

void EigenSolver::swap_vector(Vector* x, Vector* x0)
{
  Vector* swap = x;
  x  = x0;
  x0 = swap;
}

//----------------------------------------------------------------------------//
EigenSolver::size_t EigenSolver::number_iterations()
{
  return d_number_iterations;
}

} // end namespace linear_algebra

//----------------------------------------------------------------------------//
//              end of file EigenSolver.cc
//----------------------------------------------------------------------------//

