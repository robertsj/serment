//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  GlobalSolverBase.cc
 *  @brief GlobalSolverBase
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "GlobalSolverBase.hh"

namespace erme_solver
{

//----------------------------------------------------------------------------//
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

  d_local_size = d_R->number_local_rows();
  if (serment_comm::Comm::is_last()) d_local_size += 2;

  Ensure(d_residual);
}

//----------------------------------------------------------------------------//
GlobalSolverBase::~GlobalSolverBase()
{
  /* ... */
}

//----------------------------------------------------------------------------//
GlobalSolverBase::vec_dbl GlobalSolverBase::residual_norms() const
{
  return d_residual_norms;
}

//----------------------------------------------------------------------------//
void GlobalSolverBase::update_response(const double keff)
{
  /// Give the server the new keff
  d_server->update(keff);

  /// Fill the operators with updated values
  if (serment_comm::Comm::is_global())
  {
    d_R->update();
    d_F->update();
    d_A->update();
    d_L->update();
  }
}

//----------------------------------------------------------------------------//
void GlobalSolverBase::setup_initial_current(Vector &x)
{
  Require(x.local_size() >= d_R->number_local_rows());
  using serment_comm::Comm;
  if (Comm::is_global())
  {
    // Set zeroth order space/angle moments to unity, with others to zero. If
    // different energy bases are used, better starting guesses may be used.
    x.set(0.0);
    for (int i = 0; i < d_indexer->number_local_moments(); ++i)
    {
      erme_response::ResponseIndex ri = d_indexer->response_index_from_local(i);
      if (ri.azimuth + ri.polar + ri.space0 + ri.space1 == 0)
        x[i] = 1.0;
    }
    x.assemble();
    x.scale(1.0/x.norm(x.L2));
  }
}

} // end namespace erme_solver

//----------------------------------------------------------------------------//
//              end of file GlobalSolverBase.cc
//----------------------------------------------------------------------------//
