//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  GlobalSolverPicard.cc
 *  @brief GlobalSolverPicard
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "GlobalSolverPicard.hh"
#include "SteffensenUpdate.hh"
#include "AnghelUpdate.hh"
#include <cstdio>

namespace erme_solver
{

using serment_comm::Comm;

//----------------------------------------------------------------------------//
GlobalSolverPicard::GlobalSolverPicard(SP_db                db,
                                       SP_indexer           indexer,
                                       SP_server            server,
                                       SP_state             state,
                                       SP_responsecontainer responses)
  : Base(db, indexer, server, state, responses)
{

  if (Comm::is_global())
  {
    d_MR = new OperatorMR(d_R, d_M);

    // Picard-specific DB options
    int         inner_max_iters = 1e5;
    double      inner_tolerance = 1e-13;
    std::string updater         = "default";
    if (d_db->check("erme_inner_max_iters"))
      inner_max_iters = d_db->get<int>("erme_inner_max_iters");
    if (d_db->check("erme_inner_tolerance"))
      inner_tolerance = d_db->get<double>("erme_inner_tolerance");
    if (d_db->check("erme_picard_update"))
    {
      updater = d_db->get<std::string>("erme_picard_update");
    }
    if (updater == "default")
    {
      d_update = new EigenvalueUpdate();
    }
    else if (updater == "steffensen")
    {
      d_update = new SteffensenUpdate();
    }
    else if (updater == "anghel")
    {
      int anghel_scheme = AnghelUpdate::EXP;
      if (d_db->check("erme_anghel_scheme"))
        anghel_scheme = d_db->get<int>("erme_anghel_scheme");
      d_update = new AnghelUpdate(anghel_scheme);
    }
    else
    {
      THROW("Unknown Picard eigenvalue updater: " + updater);
    }

    d_innersolver =
      new linear_algebra::EigenSolver(d_MR,
                                      linear_algebra::Matrix::SP_matrix(0),
                                      inner_max_iters,
                                      inner_tolerance);

    if (serment_comm::Comm::world_rank() == 0)
      std::cout << "Using eigenvalue updater: " << updater << std::endl;

  } // end global
}

//----------------------------------------------------------------------------//
void GlobalSolverPicard::solve()
{
  Insist(serment_comm::communicator == serment_comm::world,
         "Solvers must be accessed on the world communicator.");

  if (Comm::world_rank() == 0) std::cout << " Picard: solving... " << std::endl;

  // Set the initial guess for unknowns
  SP_vector x, J;
  double    keff = 1.0;
  double    lambda = 1.0;
  if (Comm::is_global())
  {
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
  }

  // Compute initial responses
  update_response(keff);

  // Compute the initial residual norm
  double norm = 1.0;
  if (Comm::is_global())
    norm = d_residual->compute_norm(*x);
  Comm::broadcast(&norm, 1, 0);
  d_residual_norms.push_back(norm);

  // Perform outer iterations
  display(0, norm, lambda, keff);
  int it = 1;
  for (; it <= d_maximum_iterations; it++)
  {
    norm = iterate(x, keff, lambda);
    d_residual_norms.push_back(norm);
    display(it, norm, lambda, keff);
    if (norm < d_tolerance) break;
  }

  if (Comm::world_rank() == 0)
  {
    std::printf(" FINAL PICARD ITERATION EIGENVALUES: \n");
    std::printf(" **** FINAL KEFF        = %12.9f \n", keff);
    std::printf(" **** FINAL LAMBDA      = %12.9f \n", lambda);
    std::printf(" **** OUTER ITERATIONS  = %8i \n", it);
  }

  // Update the state
  d_state->update(J, keff, lambda);
}

//----------------------------------------------------------------------------//
double GlobalSolverPicard::iterate(SP_vector x, double &keff, double &lambda)
{
  if (Comm::is_global())
  {
    // Wrap current in R-sized vector
    Vector::SP_vector J(new Vector(*x, d_R->number_local_rows()));

    // Current balance (inner iterations)
    lambda = d_innersolver->solve(J);

    // Initial update of keff via gains-to-losses
    keff =  d_F->dot(*J) / (d_A->dot(*J) + d_L->leakage(*J));

    // Improved the update, if applicable
    keff = d_update->compute(keff, lambda, J);

    if (Comm::is_last)
    {
      (*x)[d_local_size - 2] = keff;
      (*x)[d_local_size - 1] = lambda;
    }
  }

  // Update responses with new keff and return the new residual
  Comm::broadcast(&keff, 1, 0);
  update_response(keff);
  double norm = 0.0;
  if (Comm::is_global()) norm = d_residual->compute_norm(*x);
  Comm::broadcast(&norm, 1, 0);
  return norm;
}

//----------------------------------------------------------------------------//
void GlobalSolverPicard::display(const size_t it,
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
//              end of file GlobalSolverPicard.cc
//----------------------------------------------------------------------------//
