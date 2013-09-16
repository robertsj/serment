//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  GlobalSolverNewton.cc
 *  @brief GlobalSolverNewton class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "GlobalSolverNewton.hh"
#include "GlobalSolverPicard.hh"
#include "FullJacobian.hh"
#include <cstdio>

namespace erme_solver
{

//----------------------------------------------------------------------------//
GlobalSolverNewton::GlobalSolverNewton(SP_db                db,
                                       SP_indexer           indexer,
                                       SP_server            server,
                                       SP_state             state,
                                       SP_responsecontainer responses)
  : Base(db, indexer, server, state, responses)
{

  d_residual = new NonlinearResidual(server, responses);

  // Delta-k for finite difference approximation.  There is a rich theory
  // on selecting this value, but something roughly equal to the square
  // root of machine epsilon typically works well.
  double eps = 1.0e-8;
  if (db->check("erme_newton_fd_epsilon"))
    eps = db->get<double>("erme_newton_fd_epsilon");

  // Optionally use the full Jacobian for everything.  This is likely less
  // efficient given the MR construction, but it might offer a bit improvement
  // in parallel
  bool full_jacobian = false;
  if (db->check("erme_newton_full_jacobian"))
    full_jacobian = 0 != db->get<int>("erme_newton_full_jacobian");

  if (full_jacobian)
    d_jacobian = new FullJacobian(d_server, d_indexer, responses, eps, false);
  else
    d_jacobian = new Jacobian(d_server, d_indexer, responses, eps);

  // Determine which matrix to use as a preconditioner.  The method can
  // be completely matrix-free in that a Jacobian is never needed explicitly.
  // However, having an explicit operator, even if just approximate, lets
  // one use factorization-based preconditioners

  std::string newton_pc = "none";
  if (db->check("erme_newton_pc"))
    newton_pc = db->get<std::string>("erme_newton_pc");

  int newton_pc_levels = 2;
  if (db->check("erme_newton_pc_levels"))
    newton_pc_levels = db->get<int>("erme_newton_pc_levels");

  int newton_pc_lag = 1;
  if (db->check("erme_newton_pc_lag"))
    newton_pc_lag = db->get<int>("erme_newton_pc_lag");

  double pc_zero_diagonal = 0.0001;
  if (db->check("erme_newton_pc_zero_diagonal"))
    pc_zero_diagonal = db->get<double>("erme_newton_pc_zero_diagonal");

  if (newton_pc == "full")
  {
    d_preconditioner =
      new FullJacobian(d_server, d_indexer, responses, eps, true, true, pc_zero_diagonal);
  }
  else if (newton_pc == "appx")
  {
    d_preconditioner =
      new FullJacobian(d_server, d_indexer, responses, eps, true, false, pc_zero_diagonal);
  }
  else
  {
    d_preconditioner = d_jacobian;
  }

  if (serment_comm::Comm::is_global())
  {
    d_solver = new linear_algebra::NonlinearSolver;
    d_solver->setup(d_db, d_residual, d_jacobian, d_preconditioner);
    if (newton_pc != "none")
    {
      d_solver->setup_pc(newton_pc_levels, newton_pc_lag);
    }
  }

}

//----------------------------------------------------------------------------//
void GlobalSolverNewton::solve()
{
  using serment_comm::Comm;
  int msg = CONTINUE;

  // Unknowns
  SP_vector x;
  double keff = 1.0;//0.9961814142140659;
  double lambda = 1.0;
  double norm = 0.0;

  // Initial guess via two picard iterations
  std::cout << "...initial guess" << std::endl;
  GlobalSolverPicard initial(d_db, d_indexer, d_server, d_state, d_responses);
  if (Comm::is_global())
  {
    Comm::set(serment_comm::global);
    x = new Vector(d_jacobian->matrix()->number_local_rows(), 0.0);
    setup_initial_current(*x);
    if (d_db->check("erme_initial_keff"))
      keff = d_db->get<double>("erme_initial_keff");
    if (Comm::is_last())
    {
      (*x)[d_local_size - 2] = keff;
      (*x)[d_local_size - 1] = lambda;
    }
    Comm::set(serment_comm::world);
  }

  std::cout << "...iterate " << std::endl;
  update_response(keff);
  //norm = d_residual->compute_norm(x.bp());
  d_residual_norms.push_back(d_residual->compute_norm(x.bp()));
  initial.iterate(x, keff, lambda);
  d_number_outer_iterations += 1;
  d_number_inner_iterations += initial.number_inner_iterations();

  std::cout << "...newton solve" << std::endl;
  if (Comm::is_global())
  {
    d_jacobian->update(x);
    d_solver->solve(x, d_tolerance, d_tolerance, d_maximum_iterations);
    msg = COMPLETED;
    Comm::broadcast(&msg, 1, 0);
    keff   = (*x)[d_local_size-2];
    lambda = (*x)[d_local_size-1];
    d_number_outer_iterations += d_solver->number_iterations();
    d_number_inner_iterations += d_solver->number_linear_iterations();
    for (int i = 0; i < d_solver->residual_norm().size(); ++i)
    {
      d_residual_norms.push_back(d_solver->residual_norm()[i]);
    }
  }
  else
  {
    while (true)
    {
      Comm::broadcast(&msg, 1, 0);
      Assert(msg == COMPLETED || msg == CONTINUE);
      if (msg == COMPLETED) break;
      // keff is broadcasted inside the call, so a dummy argument is fine
      update_response(0.0);
    }
  }

  Comm::broadcast(&keff,   1, Comm::last());
  Comm::broadcast(&lambda, 1, Comm::last());
  if (Comm::world_rank() == 0)
  {
    std::printf(" FINAL NEWTON ITERATION EIGENVALUES: \n");
    std::printf(" **** FINAL KEFF        = %12.9f \n", keff);
    std::printf(" **** FINAL LAMBDA      = %12.9f \n", lambda);
  }

  // Update the state
  d_state->update(x);
}

} // end namespace erme_solver

//----------------------------------------------------------------------------//
//              end of file GlobalSolverNewton.cc
//----------------------------------------------------------------------------//
