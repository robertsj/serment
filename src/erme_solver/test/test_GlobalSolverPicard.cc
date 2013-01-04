//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_GlobalSolverPicard.cc
 *  @brief  test_GlobalSolverPicard
 *  @author Jeremy Roberts
 *  @date   Oct 6, 2012
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_GlobalSolverPicard)

#include "erme_solver/GlobalSolverPicard.hh"
#include "utilities/TestDriver.hh"
#include "erme_response/ResponseIndexer.hh"
#include "erme_response/ResponseServer.hh"
#include "erme_geometry/NodePartitioner.hh"
#include "linear_algebra/LinearAlgebraSetup.hh"
#include <iostream>

// Setup
#include "erme_geometry/test/nodelist_fixture.hh"
using namespace erme_response;
using namespace detran_test;
using detran_utilities::soft_equiv;
using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//---------------------------------------------------------------------------//
// TEST DEFINITIONS
//---------------------------------------------------------------------------//

int test_GlobalSolverPicard(int argc, char *argv[])
{

  typedef erme_solver::GlobalSolverPicard Solver;

  //-------------------------------------------------------------------------//
  // SETUP COMM
  //-------------------------------------------------------------------------//

  typedef serment_comm::Comm Comm;

  Comm::initialize(argc, argv);
  // Setup local and global communicators.
  int number_local_comm = 1;
  if (Comm::size() > 2) number_local_comm = 2;
  Comm::setup_communicators(number_local_comm);

  //-------------------------------------------------------------------------//
  // SETUP LINEAR ALGEBRA
  //-------------------------------------------------------------------------//

  // This sets the PETSc communicator to global
  linear_algebra::initialize(argc, argv);

  //-------------------------------------------------------------------------//
  // SETUP NODES AND PARTITION
  //-------------------------------------------------------------------------//

  // Get the node list
  erme_geometry::NodeList::SP_nodelist nodes;
  if (Comm::rank() == 0)
    nodes = erme_geometry::cartesian_node_dummy_list_2d(1, 0, 0);

  // Partition the nodes
  erme_geometry::NodePartitioner partitioner;
  partitioner.partition(nodes);

  //-------------------------------------------------------------------------//
  // SETUP INDEXER AND SERVER
  //-------------------------------------------------------------------------//

  // Create parameter database
  erme_response::ResponseIndexer::SP_db db(new detran_utilities::InputDB());
  db->put<int>("dimension", 2);
  db->put<int>("erme_order_reduction", 3);

  // Create indexer
  erme_response::ResponseIndexer::SP_indexer
    indexer(new erme_response::ResponseIndexer(db, nodes));

  // Create server
  erme_response::ResponseServer::SP_server
    server(new erme_response::ResponseServer(nodes, indexer));
  server->update(1.0);

  //-------------------------------------------------------------------------//
  // SETUP STATE AND OPERATORS AND SOLVER
  //-------------------------------------------------------------------------//

  if (Comm::is_global())
  {
    // Create state
    Solver::SP_state
      state(new erme::StateERME(indexer->number_local_moments()));

    // Create response operators
    Solver::SP_M M(new erme::Connect(nodes, indexer));
    Solver::SP_R R(new erme::ResponseMatrix(nodes, indexer, server));
    Solver::SP_L L(new erme::LeakageOperator(nodes, indexer, server));
    Solver::SP_F F(new erme::FissionOperator(nodes, indexer, server));
    Solver::SP_A A(new erme::AbsorptionOperator(nodes, indexer, server));

    // Solver
    Solver::SP_solver solver(new Solver(db, indexer, server, state,
                                        R, M, F, A, L));

    solver->solve();
  }

  //-------------------------------------------------------------------------//
  // WRAP UP
  //-------------------------------------------------------------------------//

  linear_algebra::finalize();
  Comm::finalize();
  return 0;

}

//---------------------------------------------------------------------------//
//              end of file test_GlobalSolverPicard.cc
//---------------------------------------------------------------------------//
