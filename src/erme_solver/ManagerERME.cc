//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ManagerERME.cc
 *  @brief ManagerERME member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "ManagerERME.hh"
#include "comm/Comm.hh"
#include "erme_solver/GlobalSolverPicard.hh"
#include "erme_solver/GlobalSolverNewton.hh"
#include "erme_solver/GlobalSolverSteffensen.hh"
#include "linear_algebra/LinearAlgebraSetup.hh"

namespace erme_solver
{

using serment_comm::Comm;
using std::cout;
using std::endl;

//----------------------------------------------------------------------------//
ManagerERME::ManagerERME(int argc, char *argv[])
  : d_argc(argc)
  , d_argv(argv)
  , d_is_built(false)
{
  /* ... */
}

//----------------------------------------------------------------------------//
ManagerERME::SP_manager ManagerERME::Create(int argc, char *argv[])
{
  SP_manager p(new ManagerERME(argc, argv));
  return p;
}

//----------------------------------------------------------------------------//
void ManagerERME::build_comm(SP_db db, bool flag)
{
  // Note, the db is not guaranteed to be defined on all processes, as
  // only rank 0 reads from archives.
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
    if (local > Comm::size()) local = Comm::size();
  }

  // Broadcast local group count and setup comm groups
  Comm::broadcast(&local, 1, 0);
  Comm::setup_communicators(local);

  // Initialize linear algebra system
  if (flag)
    linear_algebra::initialize(d_argc, d_argv);
}
//----------------------------------------------------------------------------//
void ManagerERME::build_erme(SP_nodelist nodes)
{
  // Note, the db is not guaranteed to be defined on all processes, as
  // only rank 0 reads from archives.
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
  if (Comm::is_global())
  	d_responses = new erme::ResponseContainer(d_nodes, d_indexer, d_server);

  d_is_built = true;

  if (Comm::rank() == 0) cout << "*** ERME BUILT" << endl;
}

//----------------------------------------------------------------------------//
void ManagerERME::solve()
{
  Insist(d_is_built, "Must build the manager before using.");
  using namespace erme_solver;

  std::string solver_type = "picard";
  if (d_db->check("erme_solver_type"))
  {
    solver_type = d_db->get<std::string>("erme_solver_type");
  }

  if (solver_type == "picard")
  {
    d_solver = new GlobalSolverPicard(
      d_db, d_indexer, d_server, d_state, d_responses);
  }
  else if (solver_type == "newton")
  {
    d_solver = new GlobalSolverNewton(
      d_db, d_indexer, d_server, d_state, d_responses);
  }
  else if (solver_type == "steffensen")
  {
    d_solver = new GlobalSolverSteffensen(
      d_db, d_indexer, d_server, d_state, d_responses);
  }
  else
  {
    THROW("Unsupported solver type: " + solver_type);
  }
  d_solver->solve();
}

//----------------------------------------------------------------------------//
void ManagerERME::initialize(int argc, char *argv[])
{
  serment_comm::Comm::initialize(argc, argv);
}

//----------------------------------------------------------------------------//
void ManagerERME::finalize()
{
  linear_algebra::finalize();
  serment_comm::Comm::finalize();
}

} // end namespace erme_solver

//----------------------------------------------------------------------------//
//              end of file ManagerERME.cc
//----------------------------------------------------------------------------//
