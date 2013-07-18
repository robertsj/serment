//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Connect.cc
 *  @brief Test of Connect class
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST           \
        FUNC(test_Connect)

#include "utilities/TestDriver.hh"
#include "Connect.hh"
#include "linear_algebra/LinearAlgebraSetup.hh"
#include "erme_geometry/NodePartitioner.hh"
#include "erme_geometry/test/nodelist_fixture.hh"
#include <iostream>

// Setup

using namespace erme;
using namespace detran_test;
using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//-----------------------------------------------//
// TEST DEFINITIONS
//-----------------------------------------------//

// Test of basic public interface
int test_Connect(int argc, char *argv[])
{

  //--------------------------------------------------------------------------//
  // SETUP COMM
  //--------------------------------------------------------------------------//

  typedef serment_comm::Comm Comm;

  // Initialize comm
  Comm::initialize(argc, argv);

  // Setup local and global communicators.
  int number_local_comm = 1;
  if (Comm::size() > 2) number_local_comm = 2;
  Comm::setup_communicators(number_local_comm);

  //--------------------------------------------------------------------------//
  // SETUP LINEAR ALGEBRA
  //--------------------------------------------------------------------------//

  // This sets the PETSc communicator to global
  linear_algebra::initialize(argc, argv);

  //--------------------------------------------------------------------------//
  // SETUP NODES AND PARTITION
  //--------------------------------------------------------------------------//

  // Get the node list
  erme_geometry::NodeList::SP_nodelist nodes;
  if (Comm::rank() == 0)
    nodes = erme_geometry::cartesian_node_dummy_list_2d(1, 0, 0);

  // Partition the nodes
  erme_geometry::NodePartitioner partitioner;
  partitioner.partition(nodes);

  //--------------------------------------------------------------------------//
  // SETUP INDEXER
  //--------------------------------------------------------------------------//

  // Create parameter database
  erme_response::ResponseIndexer::SP_db db(new detran_utilities::InputDB());
  db->put<int>("dimension", 2);
  db->put<int>("erme_order_reduction", 3);

  // Create indexer
  erme_response::ResponseIndexer::SP_indexer
    indexer(new erme_response::ResponseIndexer(db, nodes));

  //--------------------------------------------------------------------------//
  // TEST
  //--------------------------------------------------------------------------//

  /*
   *  The global problem is defined only on the global processes. Hence,
   *  the operators M, R, etc. must be constructed in this communicator.
   *  Of course, the response servers live everywhere, and internally,
   *  the operators get their data from the response server root
   *  processes.
   *
   */
  if (Comm::is_global())
  {
    // Connect
    Connect M(nodes, indexer);

    // Write to binary for inspection in MATLAB
    M.display(Connect::BINARY, "connectivity.out");
  }

  //--------------------------------------------------------------------------//
  // WRAP UP
  //--------------------------------------------------------------------------//

  linear_algebra::finalize();
  Comm::finalize();
  return 0;

}

//----------------------------------------------------------------------------//
//              end of test_StateERME.cc
//----------------------------------------------------------------------------//
