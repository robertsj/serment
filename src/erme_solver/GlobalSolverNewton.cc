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
  double tol = 1.0e-8;
  if (db->check("erme_newton_fd_epsilon"))
    tol = db->get<double>("erme_newton_fd_epsilon");

  // Optionally use the full Jacobian for everything.  This is likely less
  // efficient given the MR construction, but it might offer a bit improvement
  // in parallel
  if (db->check("erme_newton_full_jacobian"))
    d_jacobian = new FullJacobian(d_server, d_indexer, responses, tol, false);
  else
    d_jacobian = new Jacobian(d_server, d_indexer, responses, tol);

  // Determine which matrix to use as a preconditioner.  The method can
  // be completely matrix-free in that a Jacobian is never needed explicitly.
  // However, having an explicit operator, even if just approximate, lets
  // one use factorization-based preconditioners
  std::string newton_pc = "none";
  if (db->check("erme_newton_pc"))
    newton_pc = db->get<std::string>("erme_newton_pc");
  if (newton_pc == "full")
  {
    d_preconditioner =
      new FullJacobian(d_server, d_indexer, responses, tol, true, true);
  }
  else if (newton_pc == "appx")
  {
    d_preconditioner =
      new FullJacobian(d_server, d_indexer, responses, tol, true, false);
  }
  else
  {
    d_preconditioner = d_jacobian;
  }

  if (serment_comm::Comm::is_global())
  {
    d_solver = new linear_algebra::NonlinearSolver;
    d_solver->setup(d_db, d_residual, d_jacobian, d_preconditioner);
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

  // Initial guess
  std::cout << "...initial guess" << std::endl;
  GlobalSolverPicard initial(d_db, d_indexer, d_server, d_state, d_responses);
  if (Comm::is_global())
  {
    if (d_db->check("erme_initial_keff"))
      keff = d_db->get<double>("erme_initial_keff");
    x = new Vector(d_jacobian->matrix()->number_local_rows(), 0.0);
    setup_initial_current(*x);
  }
  std::cout << "...iterate " << std::endl;
  update_response(keff);
  norm = initial.iterate(x, keff, lambda);
  norm = initial.iterate(x, keff, lambda);
  display_response("jac");
  x->display(x->BINARY, "X.out");
//  update_response(keff);
//  norm = initial.iterate(x, keff, lambda);
//  update_response(keff);
//  norm = initial.iterate(x, keff, lambda);
//  update_response(keff);
//  norm = initial.iterate(x, keff, lambda);

  // Solve
  std::cout << "...newton solve" << std::endl;
  if (Comm::is_global())
  {
    d_jacobian->update(x);
    d_jacobian->matrix()->display(d_jacobian->matrix()->BINARY, "JAC.OUT");
    double nx = x->norm();
    //d_preconditioner->update(x);
    Assert(nx == x->norm());
    //d_preconditioner->matrix()->display(d_jacobian->matrix()->BINARY, "PRE.OUT");
    d_solver->solve(x);
    msg = COMPLETED;
    Comm::broadcast(&msg, 1, 0);
    keff   = (*x)[d_local_size-2];
    lambda = (*x)[d_local_size-1];
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
  x->display(x->BINARY, "Xfinal.out");

  Comm::broadcast(&keff,   1, Comm::last());
  Comm::broadcast(&lambda, 1, Comm::last());
  if (Comm::world_rank() == 0)
  {
    std::printf(" FINAL NEWTON ITERATION EIGENVALUES: \n");
    std::printf(" **** FINAL KEFF        = %12.9f \n", keff);
    std::printf(" **** FINAL LAMBDA      = %12.9f \n", lambda);
  }

  // Update the state
  d_state->set_k(keff);
  //d_state->update(J, keff, lambda);
}

} // end namespace erme_solver

//----------------------------------------------------------------------------//
//              end of file GlobalSolverNewton.cc
//----------------------------------------------------------------------------//
