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

#include "TestDriver.hh"
#include "ResponseMatrix.hh"
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

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

// Test of basic public interface
int test_ResponseMatrix(int argc, char *argv[])
{
  typedef serment_comm::Comm Comm;

  Comm::initialize(argc, argv);
  PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

  // Get the node list
  erme_geometry::NodeList::SP_nodelist nodes;
  if (Comm::rank() == 0)
    nodes = erme_geometry::cartesian_node_dummy_list_2d();

  // Partition the nodes
  erme_geometry::NodePartitioner partitioner;
  partitioner.partition(nodes);

  // Create parameter database
  erme_response::ResponseIndexer::SP_db db(new detran::InputDB());
  db->put<int>("dimension", 2);
  db->put<int>("erme_order_reduction", 3);

  // Create indexer
  erme_response::ResponseIndexer::SP_indexer
    indexer(new erme_response::ResponseIndexer(db, nodes));

  // Creater server
  erme_response::ResponseServer::SP_server
    server(new erme_response::ResponseServer(nodes, indexer));

  // Create new scope for M, since M must be destructed before finalization.
  {
    // ResponseMatrix
    ResponseMatrix R(nodes, indexer, server);

    // Update.  Since this uses the dummy data, this simply fills R.
    R.update(1.0);

    // Write to binary for inspection in MATLAB
    R.display(ResponseMatrix::BINARY, "response_matrix.out");
  }

  PetscFinalize();
  Comm::finalize();
  return 0;

}

//---------------------------------------------------------------------------//
//              end of test_StateERME.cc
//---------------------------------------------------------------------------//
