//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  GlobalSolverNewton.cc
 *  @brief GlobalSolverNewton class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "GlobalSolverNewton.hh"

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

  d_residual = new NonlinearResidual(d_R, d_M, d_F, d_A, d_L);
  d_jacobian = new Jacobian(d_server, responses, 1.0e-8);

  d_solver = new linear_algebra::NonlinearSolver;
  d_solver->setup(d_db, d_residual, d_jacobian, d_jacobian);
}

//----------------------------------------------------------------------------//
void GlobalSolverNewton::solve()
{
  SP_vector x(new Vector(d_jacobian->matrix()->number_local_rows(), 0.0));

  setup_initial_current(*x);

  d_solver->solve(x);
}

} // end namespace erme_solver
