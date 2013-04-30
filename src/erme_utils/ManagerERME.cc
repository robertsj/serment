//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   ManagerERME.cc
 *  @brief  ManagerERME member definitions
 *  @author Jeremy Roberts
 *  @date   Oct 4, 2012
 */
//----------------------------------------------------------------------------//

#include "ManagerERME.hh"
#include "comm/Comm.hh"
#include "erme_solver/GlobalSolverPicard.hh"
#include "linear_algebra/LinearAlgebraSetup.hh"

namespace erme_utils
{

//----------------------------------------------------------------------------//
ManagerERME::ManagerERME(int argc, char *argv[])
  : d_argc(argc)
  , d_argv(argv)
  , d_is_built(false)
{
  serment_comm::Comm::initialize(argc, argv);
}

//----------------------------------------------------------------------------//
ManagerERME::ManagerERME(int argc, char *argv[], SP_db db)
  : d_argc(argc)
  , d_argv(argv)
  , d_db(db)
  , d_is_built(false)
{
  serment_comm::Comm::initialize(d_argc, d_argv);
  build_comm(db);
}

//----------------------------------------------------------------------------//
ManagerERME::SP_manager ManagerERME::Create(int argc, char *argv[], SP_db db)
{
  SP_manager p(new ManagerERME(argc, argv, db));
  return p;
}

//----------------------------------------------------------------------------//
void ManagerERME::build_comm(SP_db db)
{
  // Note, the db is not guaranteed to be defined on all processes, as
  // only rank 0 reads from archives.

  using serment_comm::Comm;
  using std::cout;
  using std::endl;

  if (Comm::rank() == 0)
  {
    Insist(db, "Parameter database is NULL");
  }
  d_db = db;
  Comm::broadcast(d_db, 0);
  Assert(d_db);

  if (Comm::rank() == 0) cout << "BUILDING MANAGER" << endl;

  // Check the desired local group structure
  int local = 1;
  if (Comm::world_rank() == 0)
  {
    if (d_db->check("comm_local_groups"))
    {
      local = d_db->get<int>("comm_local_groups");
    }
    Assert(local > 0);
  }

  // Broadcast local group count and setup comm groups
  Comm::broadcast(&local, 1, 0);
  Comm::setup_communicators(local);

  // Initialize linear algebra system
  linear_algebra::initialize(d_argc, d_argv);
}
//----------------------------------------------------------------------------//
void ManagerERME::build_erme(SP_nodelist nodes)
{
  // Note, the db is not guaranteed to be defined on all processes, as
  // only rank 0 reads from archives.

  using std::cout;
  using std::endl;
  using serment_comm::Comm;

  if (Comm::rank() == 0)
  {
    Insist(Comm::is_comm_built(),
           "Communicators must be built before problem construction.")
    Insist(nodes, "Node list is NULL");
  }
  d_nodes = nodes;

  if (Comm::rank() == 0) cout << "*** BUILDING ERME " << endl;

  // Check if this is a database problem
  std::string dbname = "";
  size_t dborder     = 1;
  if (Comm::rank() == 0)
  {
    if (d_db->check("response_db_name"))
    {
      dbname = d_db->get<std::string>("response_db_name");
    }
    if (d_db->check("response_db_order"))
    {
      dborder = d_db->get<int>("response_db_order");
    }
  }

  // Create partitioner and partition
  Partitioner partitioner;
  partitioner.partition(d_nodes);

  // Create indexer
  if (Comm::rank() == 0) cout << "****** BUILDING INDEXER" << endl;
  d_indexer = new erme_response::ResponseIndexer(d_db, d_nodes);

  // Create server
  if (Comm::rank() == 0) cout << "****** BUILDING SERVER" << endl;
  d_server = new erme_response::
    ResponseServer(d_nodes, d_indexer, dbname, dborder);

  // Create state
  if (Comm::rank() == 0) cout << "****** BUILDING STATE" << endl;
  d_state = new erme::StateERME(d_indexer->number_local_moments());

  // Create response operators
  if (Comm::rank() == 0) cout << "****** BUILDING OPERATORS" << endl;

  // Operators and live on global communicator
  if (Comm::is_global())
  	d_responses = new erme::ResponseContainer(d_nodes, d_indexer, d_server);
  d_is_built = true;

}

//----------------------------------------------------------------------------//
void ManagerERME::solve()
{
  // Preconditions
  Insist(d_is_built, "Must build the manager before using.");

  // Create solver
  d_solver = new erme_solver::GlobalSolverPicard(
    d_db, d_indexer, d_server, d_state, d_responses);

  // Solve
  d_solver->solve();
}

//----------------------------------------------------------------------------//
void ManagerERME::finalize()
{
  /* ... */
}

} // end namespace erme_utils

//----------------------------------------------------------------------------//
//              end of file ManagerERME.cc
//----------------------------------------------------------------------------//
