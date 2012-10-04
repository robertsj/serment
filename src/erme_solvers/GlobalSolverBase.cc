//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   GlobalSolverBase.cc
 *  @brief  GlobalSolverBase
 *  @author Jeremy Roberts
 *  @date   Oct 1, 2012
 */
//---------------------------------------------------------------------------//

#include "GlobalSolverBase.hh"

namespace erme_solver
{

GlobalSolverBase::
GlobalSolverBase(SP_db db, SP_server server, SP_state state,
                 SP_R R, SP_M M, SP_F F, SP_A A, SP_L L)
  : d_db(db)
  , d_server(server)
  , d_state(state)
  , d_R(R)
  , d_M(M)
  , d_F(F)
  , d_A(A)
  , d_L(L)
{
  // Preconditions
  Require(d_db);
  Require(d_server);
  Require(d_state);
  Require(d_R);
  Require(d_M);
  Require(d_F);
  Require(d_A);
  Require(d_L);

  // Create residual
  d_residual = new NonlinearResidual(R, M, F, A, L);

  // Set convergence criteria


  // Postconditions
  Ensure(d_residual);
}

GlobalSolverBase::~GlobalSolverBase()
{ /* ... */ }

} // end namespace erme_solver

//---------------------------------------------------------------------------//
//              end of file GlobalSolverBase.cc
//---------------------------------------------------------------------------//
