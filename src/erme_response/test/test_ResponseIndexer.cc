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

  // no order reduction
  {
    db->put<int>("erme_order_reduction", 0);
    ResponseIndexer indexer(db, nodes);

    TEST(indexer.number_nodes()               == 4);
    TEST(indexer.number_node_moments(0)       == 180);
    TEST(indexer.node_index(0, 0, 0).azimuth  == 0);
    TEST(indexer.node_index(0, 3, 44).polar   == 2);
    TEST(indexer.node_index(0, 3, 44).azimuth == 2);
    TEST(indexer.node_index(0, 3, 44).space0  == 4);
  }

  // space-angle reduction
  {

    db->put<int>("erme_order_reduction", 3);
    ResponseIndexer indexer(db, nodes);
    TEST(indexer.number_nodes()                 == 4);
    TEST(indexer.number_node_moments(0)         == 108);
    TEST(indexer.node_index(0, 0, 0).azimuth    == 0);
    TEST(indexer.node_index(0, 0, 26).polar     == 2);
    TEST(indexer.node_index(0, 0, 26).azimuth   == 2);
    TEST(indexer.node_index(0, 0, 26).space0    == 0);

    // Nodes
    for (int n = nodes.lower_bound(); n < nodes.upper_bound(); n++)
    {
      // Moment index on node
      int n_index = 0;
      // Surfaces
      for (int s = 0; s < 4; s++)
      {
        // Neighbor and its surface
        int neigh_n = nodes.neighbor(n, s).neighbor();
        int neigh_s = nodes.neighbor(n, s).surface();
        // All moments on surface s
        for (int m = 0; m < indexer.number_surface_moments(n, s); m++, n_index++)
        {
          if (nodes.neighbor(n, s).neighbor() >= 0)
          {
            ResponseIndex neigh_index = indexer.node_index(neigh_n, neigh_s, m);
            ResponseIndex index       = indexer.node_index(n, s, m);
            TEST(neigh_index.surface  == neigh_s);
            TEST(neigh_index.energy   == index.energy);
            TEST(neigh_index.polar    == index.polar);
            TEST(neigh_index.azimuth  == index.azimuth);
            TEST(neigh_index.space0   == index.space0);
            TEST(neigh_index.space1   == index.space1);
          }
          else
          {
            TEST(nodes.neighbor(n, s).neighbor() == erme_geometry::Node::VACUUM);
          }
        } // end moments
      } // surfaces
    } // nodes

  }

  Comm::finalize();
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_StateERME.cc
//---------------------------------------------------------------------------//
