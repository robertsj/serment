//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  GlobalSolverSteffensen.cc
 *  @brief GlobalSolverSteffensen member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "GlobalSolverSteffensen.hh"
#include "GlobalSolverPicard.hh"
#include <cstdio>

namespace erme_solver
{

using serment_comm::Comm;

//----------------------------------------------------------------------------//
GlobalSolverSteffensen::GlobalSolverSteffensen(SP_db                db,
                                               SP_indexer           indexer,
                                               SP_server            server,
                                               SP_state             state,
                                               SP_responsecontainer responses)
  : Base(db, indexer, server, state, responses)
{

}

//----------------------------------------------------------------------------//
void GlobalSolverSteffensen::solve()
{
  using serment_comm::Comm;
  int msg = CONTINUE;

  // Unknowns
  SP_vector x, J;
  double keff   = 1.0; //0.9961814142140659;
  double lambda = 1.0;
  double norm   = 0.0;

  // Picard solver
  GlobalSolverPicard picard(d_db, d_indexer, d_server, d_state, d_responses);

  // Initial guess
  std::cout << "...initial guess" << std::endl;
  if (Comm::is_global())
  {
    Comm::set(serment_comm::global);
    x = new Vector(d_local_size, 1.0);
    setup_initial_current(*x);
    if (d_db->check("erme_initial_keff"))
      keff = d_db->get<double>("erme_initial_keff");
    if (Comm::is_last())
    {
      (*x)[d_local_size - 2] = keff;
      (*x)[d_local_size - 1] = lambda;
    }
    J = new Vector(*x, d_R->number_local_rows());
    Comm::set(serment_comm::world);
  }
  int count = 0;
  norm = d_residual->compute_norm(x.bp()); ++count;
  d_residual_norms.push_back(norm);

  // Perform iterations
  display(0, norm, lambda, keff);
  int it = 1;
  for (; it <= d_maximum_iterations; it++)
  {
    double k0 = keff, k1, k2;
    picard.iterate(x, k1, lambda); ++count;
    picard.iterate(x, k2, lambda); ++count;
    keff = k0 - (k0 - k1) * (k0 - k1) / (k2 - 2.0 * k1 + k0);
    norm = d_residual->compute_norm(x.bp());
    d_residual_norms.push_back(norm);
    display(it, norm, lambda, keff);
    if (norm < d_tolerance) break;
  }
  std::cout << " COUNT = " << count << std::endl;
  // Update the state
  d_state->update(J, keff, lambda);
}

//----------------------------------------------------------------------------//
void GlobalSolverSteffensen::display(const size_t it,
                                     const double norm,
                                     const double keff,
                                     const double lambda)
{
  if (serment_comm::Comm::world_rank() == 0)
  {
    printf(" %3i %12.8e %20.16f %20.16f \n", it, norm, lambda, keff);
  }
}


} // end namespace erme_solver

//----------------------------------------------------------------------------//
//              end of file GlobalSolverSteffensen.cc
//----------------------------------------------------------------------------//
