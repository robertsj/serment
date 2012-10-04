//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   GlobalSolverPicard.i.hh
 *  @brief  GlobalSolverPicard inline member definitions
 *  @author Jeremy Roberts
 *  @date   Oct 1, 2012
 */
//---------------------------------------------------------------------------//

#ifndef erme_solver_GLOBALSOLVERPICARD_I_HH_
#define erme_solver_GLOBALSOLVERPICARD_I_HH_

#include <cstdio>

namespace erme_solver
{

void GlobalSolverPicard::solve()
{

  // Hard-coded initial guess
  double keff   = 1.0;
  double lambda = 1.0;

  // Initialize balance parameters
  double loss       = 0;
  double gain       = 0;
  double absorption = 0;
  double leakage    = 0;;

  // Count total iterations
  int innertot = 0;

  // Compute the initial residual norm
  double norm = d_residual->compute_norm(*d_J0, keff, lambda);

  //-------------------------------------------------------------------------//
  // OUTER ITERATIONS
  //-------------------------------------------------------------------------//

  int it = 1; // count outer iteration
  for (; it <= d_maximum_iterations; it++)
  {

    //-----------------------------------------------------------------------//
    // INNER ITERATIONS -- solves M*R*X = lambda*X
    //-----------------------------------------------------------------------//

    lambda = d_innersolver->solve(d_J0);

    //-----------------------------------------------------------------------//
    // EIGENVALUE UPDATE -- k = fission / (absorption + leakage)
    //-----------------------------------------------------------------------//

    // Compute gains and losses.
    gain        = d_F->dot(*d_J0);
    absorption  = d_A->dot(*d_J0);
    leakage     = d_L->leakage(*d_J0);
    loss        = absorption + leakage;

    // Update keff
    keff = d_update->compute(keff, d_J0);

    // Update the responses
    update_response(keff);

    // Compute the norm of the nonlinear residual.
    norm = d_residual->compute_norm(*d_J0, keff, lambda);

    printf(" PICARD IT: %3i NORM: %8.6e LAMBDA: %12.9f KEFF: %12.9f INNERS %8i \n",
           it, norm, lambda, keff, innertot);

    // Check convergence
    if (norm < d_tolerance) break;

  } // end outers

  std::printf(" FINAL POWER ITERATION EIGENVALUES: \n");
  std::printf(" **** FINAL KEFF        = %12.9f \n", keff);
  std::printf(" **** FINAL LAMBDA      = %12.9f \n", lambda);
  std::printf(" **** OUTER ITERATIONS  = %8i \n", it);
  std::printf(" **** INNER ITERATIONS  = %8i \n", innertot);

  // Update the state

  return;

}

} // end namespace erme_solver

#endif // erme_solver_GLOBALSOLVERPICARD_I_HH_

//---------------------------------------------------------------------------//
//              end of file GlobalSolverPicard.i.hh
//---------------------------------------------------------------------------//
