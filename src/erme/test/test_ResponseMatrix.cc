//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_ResponseMatrix.cc
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 * \brief  Test of ResponseMatrix class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST           \
        FUNC(test_ResponseMatrix)

#include "utilities/TestDriver.hh"
#include "ResponseMatrix.hh"
#include "linear_algebra/LinearAlgebraSetup.hh"
#include "erme_geometry/NodePartitioner.hh"
#include "erme_geometry/test/nodelist_fixture.hh"
#include <iostream>

// Setup

using namespace erme;
using namespace detran_test;
using detran_utilities::soft_equiv;
using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

// Test of basic public interface
int test_ResponseMatrix(int argc, char *argv[])
{

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

  // Creater server
  erme_response::ResponseServer::SP_server
    server(new erme_response::ResponseServer(nodes, indexer));

  /*
   *  Update the server, which leads to the first local computations.  In
   *  some cases, this might be the only true computation.  For instance,
   *  if we define the responses as an expansion (R = R0 + R1/k + ...), then
   *  this would compute the expansion coefficients once, store them
   *  locally, and interpolate/expand on-the-fly.
   */
  server->update(1.0);

  // Create R only on the global communicator
  if (Comm::is_global())
  {
    // ResponseMatrix
    ResponseMatrix R(nodes, indexer, server);

    // Update.  Since this uses the dummy data, this simply fills R.
    R.update();

    // Write to binary for inspection in MATLAB
    R.display(ResponseMatrix::BINARY, "response_matrix.out");

  }

  //-------------------------------------------------------------------------//
  // WRAP UP
  //-------------------------------------------------------------------------//

  linear_algebra::finalize();
  Comm::finalize();
  return 0;

}

//---------------------------------------------------------------------------//
//              end of test_StateERME.cc
//---------------------------------------------------------------------------//
