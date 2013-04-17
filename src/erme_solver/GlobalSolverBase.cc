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

//---------------------------------------------------------------------------//
GlobalSolverBase::
GlobalSolverBase(SP_db 									db,
                 SP_indexer 						indexer,
                 SP_server 							server,
                 SP_state 							state,
                 SP_responsecontainer 	responses)
  : d_db(db)
  , d_indexer(indexer)
  , d_server(server)
  , d_state(state)
  , d_maximum_iterations(100)
  , d_tolerance(1.0e-6)
{
	using serment_comm::Comm;

  // Preconditions
  Require(d_db);
  Require(d_indexer);
  Require(d_server);
  Require(d_state);
  if (serment_comm::Comm::is_global())
  {
  	Require(responses);
  	d_M = responses->M;
  	d_R = responses->R;
  	d_L = responses->L;
  	d_A = responses->A;
  	d_F = responses->F;
  	Ensure(d_M);
  	Ensure(d_R);
  	Ensure(d_L);
  	Ensure(d_A);
  	Ensure(d_F);
  }

  // Create residual
  d_residual = new NonlinearResidual(d_R, d_M, d_F, d_A, d_L);

  // Set convergence criteria
  if (d_db->check("erme_maximum_iterations"))
  {
    d_maximum_iterations = d_db->get<int>("erme_maximum_iterations");
    Assert(d_maximum_iterations >= 0);
  }
  if (d_db->check("erme_tolerance"))
  {
    d_tolerance = d_db->get<double>("erme_tolerance");
    Assert(d_tolerance >= 0.0);
  }

  // Postconditions
  Ensure(d_residual);
}

//---------------------------------------------------------------------------//
GlobalSolverBase::~GlobalSolverBase()
{ /* ... */ }

//---------------------------------------------------------------------------//
GlobalSolverBase::vec_dbl GlobalSolverBase::residual_norms() const
{
  return d_residual_norms;
}

} // end namespace erme_solver

//---------------------------------------------------------------------------//
//              end of file GlobalSolverBase.cc
//---------------------------------------------------------------------------//
