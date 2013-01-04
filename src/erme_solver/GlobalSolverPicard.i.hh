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

//---------------------------------------------------------------------------//
void GlobalSolverPicard::solve()
{
  std::cout << " Picard: solving... " << std::endl;

  // Get initial keff guess
  double keff = 1.0;
  if (d_db->check("erme_initial_keff"))
    keff = d_db->get<double>("erme_initial_keff");
  // Start lambda at unity
  double lambda = 1.0;
  // Set space-angle zeroth order moments to uniform guess.  Note
  // this doesn't actually get used by SLEPc but would will in a
  // hand-coded power iteration.
  for (int i = 0; i < d_indexer->number_local_moments(); ++i)
  {
    erme_response::ResponseIndex ri = d_indexer->response_index_from_local(i);
    if (ri.azimuth + ri.polar + ri.space0 + ri.space1 == 0) (*d_J0)[i] = 1.0;
  }
  // Ensure a normalized initial guess and initialize the responses
  d_J0->scale(1.0/d_J0->norm(d_J0->L2));

  update_response(keff);

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

//    std::cout << " GAIN = "     << gain << std::endl;
//    std::cout << " ABS = "      << absorption << std::endl;
//    std::cout << " leakage = "  << leakage << std::endl;

//    d_R->display(d_L->BINARY, "R.out");
//    d_M->display(d_L->BINARY, "M.out");
//    d_L->display(d_L->BINARY, "L.out");
//    d_L->display_leakage();

    // Initial update of keff
    keff = gain / loss;

    // Improved the update, if applicable
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
