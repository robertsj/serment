//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ManagerERME.cc
 *  @brief  ManagerERME member definitions
 *  @author Jeremy Roberts
 *  @date   Oct 4, 2012
 */
//---------------------------------------------------------------------------//

#include "ManagerERME.hh"
#include "comm/Comm.hh"
#include "erme_solver/GlobalSolverPicard.hh"
#include "linear_algebra/LinearAlgebraSetup.hh"

namespace erme_utils
{

//---------------------------------------------------------------------------//
ManagerERME::ManagerERME(int argc, char *argv[], SP_db db)
  : d_db(db)
  , d_is_built(false)
{
  // Preconditions
  Require(db);

  std::cout << "Constructing ManagerERME" << std::endl;

  // Initialize comm
  serment_comm::Comm::initialize(argc, argv);

  // Check the desired local group structure
  int local = 1;
  if (serment_comm::Comm::rank() == 0)
  {
    if (d_db->check("comm_local_groups"))
    {
      local = d_db->get<int>("comm_local_groups");
    }
    Assert(local > 0);
  }

  // Broadcast local group count and setup comm groups
  serment_comm::Comm::broadcast(&local, 1, 0);
  serment_comm::Comm::setup_communicators(local);

  // Initialize linear algebra system
  linear_algebra::initialize(argc, argv);

}

//---------------------------------------------------------------------------//
ManagerERME::SP_manager ManagerERME::Create(int argc, char *argv[], SP_db db)
{
  SP_manager p(new ManagerERME(argc, argv, db));
  return p;
}

//---------------------------------------------------------------------------//
void ManagerERME::build_erme(SP_nodelist nodes)
{
  // Preconditions
  Require(nodes);
  d_nodes = nodes;

  using std::cout;
  using std::endl;

  cout << "*** BUILDING ERME " << endl;

  // Check if this is a database problem
  std::string dbname = "";
  size_t dborder = 1;
  if (d_db->check("response_db_name"))
  {
    dbname = d_db->get<std::string>("response_db_name");
  }
  if (d_db->check("response_db_order"))
  {
    dborder = d_db->get<int>("response_db_order");
  }

  // Create partitioner and partition
  Partitioner partitioner;
  partitioner.partition(d_nodes);

  // Create indexer
  cout << "****** BUILDING INDEXER" << endl;
  d_indexer = new erme_response::ResponseIndexer(d_db, d_nodes);

  // Create server
  cout << "****** BUILDING SERVER" << endl;
  d_server = new erme_response::
    ResponseServer(d_nodes, d_indexer, dbname, dborder);

  // Create state
  cout << "****** BUILDING STATE" << endl;
  d_state = new erme::StateERME(d_indexer->number_local_moments());

  // Create response operators
  cout << "****** BUILDING OPERATORS" << endl;
  d_M = new erme::Connect(d_nodes, d_indexer);
  d_R = new erme::ResponseMatrix(d_nodes, d_indexer, d_server);
  d_L = new erme::LeakageOperator(d_nodes, d_indexer, d_server);
  d_F = new erme::FissionOperator(d_nodes, d_indexer, d_server);
  d_A = new erme::AbsorptionOperator(d_nodes, d_indexer, d_server);
  d_is_built = true;
}

//---------------------------------------------------------------------------//
void ManagerERME::solve()
{
  // Preconditions
  Insist(d_is_built, "Must build the manager before using.");

  // Create solver
  d_solver = new erme_solver::GlobalSolverPicard(
    d_db, d_indexer, d_server, d_state, d_R, d_M, d_F, d_A, d_L);

  // Solve
  d_solver->solve();
}

//---------------------------------------------------------------------------//
void ManagerERME::finalize()
{
  /* ... */
}

} // end namespace erme_utils

//---------------------------------------------------------------------------//
//              end of file ManagerERME.cc
//---------------------------------------------------------------------------//
