//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_ResponseIndexer.cc
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 * \brief  Test of ResponseIndexer class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST           \
        FUNC(test_ResponseIndexer)

#include "TestDriver.hh"
#include "erme_response/ResponseIndexer.hh"
#include "erme_geometry/NodePartitioner.hh"
#include <iostream>

// Setup
#include "erme_geometry/test/nodelist_fixture.hh"

using namespace erme_response;
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
int test_ResponseIndexer(int argc, char *argv[])
{
  typedef serment_comm::Comm Comm;

  Comm::initialize(argc, argv);

  // Get the node list
  erme_geometry::NodeList nodes;
  if (Comm::rank() == 0)
    nodes = erme_geometry::cartesian_node_detran_list_2d();

  // Partition the nodes
  erme_geometry::NodePartitioner partitioner;
  partitioner.partition(nodes);

  // Create parameter database
  ResponseIndexer::SP_db db(new detran::InputDB());
  db->put<int>("dimension", 2);

  int num_nodes = 1;
  if (Comm::size() == 1)
  {
    num_nodes = 4;
  }
  if (Comm::size() == 2)
  {
    num_nodes = 2;
  }
  if (Comm::size() == 3)
  {
    if (Comm::rank() == 2) num_nodes = 2;
  }
  if (Comm::size() == 4)
  {
    if (Comm::rank() < 4)
      num_nodes = 1;
    else
      num_nodes = 0;
  }
  return 0;
  // no order reduction
  {
    db->put<int>("erme_order_reduction", 0);
    ResponseIndexer indexer(db, nodes);

    TEST(indexer.number_nodes()         == num_nodes);
    TEST(indexer.number_moments(0)      == 180);
    TEST(indexer.index(0, 0).azimuth    == 0);
    TEST(indexer.index(179, 0).polar    == 2);
    TEST(indexer.index(179, 0).azimuth  == 2);
    TEST(indexer.index(179, 0).space0   == 4);
  }

  // space-angle reduction
  {
    db->put<int>("erme_order_reduction", 3);
    ResponseIndexer indexer(db, nodes);
    TEST(indexer.number_nodes()         == num_nodes);
    TEST(indexer.number_moments(0)      == 108);
    TEST(indexer.index(0, 0).azimuth    == 0);
    TEST(indexer.index(107, 0).polar    == 2);
    TEST(indexer.index(107, 0).azimuth  == 2);
    TEST(indexer.index(107, 0).space0   == 0);
  }

  Comm::finalize();
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_StateERME.cc
//---------------------------------------------------------------------------//
