//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Manager.cc
 *  @brief  Manager member definitions
 *  @author Jeremy Roberts
 *  @date   Oct 4, 2012
 */
//---------------------------------------------------------------------------//

#include "Manager.hh"
//
#include "erme_solver/GlobalSolverPicard.hh"

namespace erme_utils
{

Manager::Manager(SP_db db)
  : d_db(db)
{
  // Preconditions
  Require(db);
}

Manager::~Manager(){}

Manager::SP_manager Manager::Create(SP_db db)
{
  SP_manager p(new Manager(db));
  return p;
}

void Manager::build_erme(SP_nodelist nodes)
{
  // Preconditions
  Require(nodes);
  d_nodes = nodes;

  // Check the desired local group structure
  int local = 1;
  if (d_db->check("comm_local_groups"))
  {
    local = d_db->get<int>("comm_local_groups");
  }

  // Create partitioner and partition
  Partitioner partitioner;
  partitioner.partition(d_nodes);

  // Create indexer
  d_indexer = new erme_response::ResponseIndexer(d_db, d_nodes);

  // Create server
  d_server = new erme_response::ResponseServer(d_nodes, d_indexer);

  // Create state
  d_state = new erme::StateERME(d_indexer->number_local_moments());

  // Create response operators
  d_M = new erme::Connect(d_nodes, d_indexer);
  d_R = new erme::ResponseMatrix(d_nodes, d_indexer, d_server);
  d_L = new erme::LeakageOperator(d_nodes, d_indexer, d_server);
  d_F = new erme::FissionOperator(d_nodes, d_indexer, d_server);
  d_A = new erme::AbsorptionOperator(d_nodes, d_indexer, d_server);

}

void Manager::solve()
{
  // Create solver
  d_solver = new erme_solver::GlobalSolverPicard(
    d_db, d_server, d_state, d_R, d_M, d_F, d_A, d_L);

  // Solve
  d_solver->solve();

}



} // end namespace erme_utils

//---------------------------------------------------------------------------//
//              end of file Manager.cc
//---------------------------------------------------------------------------//
