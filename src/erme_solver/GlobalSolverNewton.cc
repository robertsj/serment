//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   GlobalSolverNewton.cc
 *  @author robertsj
 *  @date   Dec 6, 2012
 *  @brief  GlobalSolverNewton class definition.
 */
//---------------------------------------------------------------------------//

#include "GlobalSolverNewton.hh"

namespace erme_solver
{

GlobalSolverNewton::GlobalSolverNewton(SP_db 			db,
		                                   SP_indexer indexer,
		                                   SP_server 	server,
		                                   SP_state 	state,
		                                   SP_R 			R,
		                                   SP_M 			M,
		                                   SP_F 			F,
		                                   SP_A 			A,
		                                   SP_L 			L)
  : Base(db, indexer, server, state, R, M, F, A, L)
{

}

void GlobalSolverNewton::solve()
{

}

} // end namespace erme_solver
