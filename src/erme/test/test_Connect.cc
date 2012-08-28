//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_Connect.cc
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 * \brief  Test of Connect class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST           \
        FUNC(test_Connect)

#include "TestDriver.hh"
#include "Connect.hh"
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
  PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
  RUN(argc, argv);
  PetscFinalize();
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

// Test of basic public interface
int test_Connect(int argc, char *argv[])
{
  typedef serment_comm::Comm Comm;

  // Get the node list
  erme_geometry::NodeList nodes;
  if (Comm::rank() == 0)
    nodes = erme_geometry::cartesian_node_detran_list_2d();

  // Partition the nodes
  erme_geometry::NodePartitioner partitioner;
  partitioner.partition(nodes);

  // Create parameter database
  erme_response::ResponseIndexer::SP_db db(new detran::InputDB());
  db->put<int>("dimension", 2);

  db->put<int>("erme_order_reduction", 0);
  erme_response::ResponseIndexer indexer(db, nodes);

  // Connect
  Connect(nodes, indexer);

  return 0;

}

//---------------------------------------------------------------------------//
//              end of test_StateERME.cc
//---------------------------------------------------------------------------//
