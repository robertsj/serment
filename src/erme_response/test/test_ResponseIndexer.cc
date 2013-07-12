//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   test_ResponseIndexer.cc
 *  @author Jeremy Roberts
 *  @date   Aug 19, 2012
 *  @brief  Test of ResponseIndexer class.
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                               \
        FUNC(test_ResponseIndexer)              \
        FUNC(test_ResponseIndexer_zeroth_order)

#include "utilities/TestDriver.hh"
#include "erme_response/ResponseIndexer.hh"
#include "erme_geometry/NodePartitioner.hh"
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

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

// Test of basic public interface
int test_ResponseIndexer(int argc, char *argv[])
{

  //--------------------------------------------------------------------------//
  // SETUP COMM
  //--------------------------------------------------------------------------//

  typedef serment_comm::Comm Comm;

  Comm::initialize(argc, argv);
  if (Comm::size() < 3)
    Comm::setup_communicators(1);
  else
    Comm::setup_communicators(2);

  //--------------------------------------------------------------------------//
  // CREATE NODE LIST AND PARTITION
  //--------------------------------------------------------------------------//

  // Get the node list
  erme_geometry::NodeList::SP_nodelist nodes;
  if (Comm::rank() == 0)
    nodes = erme_geometry::cartesian_node_dummy_list_2d(4, 2, 2);

  // Partition the nodes
  erme_geometry::NodePartitioner partitioner;
  partitioner.partition(nodes);

  // Create parameter database
  ResponseIndexer::SP_db db(new detran_utilities::InputDB());
  db->put<int>("dimension", 2);

  //--------------------------------------------------------------------------//
  // TEST INDEXER WITH NO ORDER REDUCTION
  //--------------------------------------------------------------------------//

  {
    db->put<int>("erme_order_reduction", 0);
    ResponseIndexer indexer(db, nodes);

    TEST(indexer.number_nodes()                   == 4);
    TEST(indexer.number_node_moments(0)           == 180);
    TEST(indexer.response_index(0, 0, 0).azimuth  == 0);
    TEST(indexer.response_index(0, 3, 44).polar   == 2);
    TEST(indexer.response_index(0, 3, 44).azimuth == 2);
    TEST(indexer.response_index(0, 3, 44).space0  == 4);

    if (Comm::size() < 3)
    {
      // One local group
      if (Comm::rank() == 0)
      {
        TEST(indexer.local_index_to_global(  0) ==   0);
        TEST(indexer.local_index_to_global(180) == 180);
        TEST(indexer.local_index_to_global(360) == 360);

      }
    }
    else
    {
      // Two local groups
      if (Comm::rank() == 0 and Comm::local_group() == 0)
      {
        TEST(indexer.local_index_to_global(  0) ==   0);
        TEST(indexer.local_index_to_global(180) == 180);
      }
      if (Comm::rank() == 0 and Comm::local_group() == 1)
      {
        TEST(indexer.local_index_to_global(  0) == 360);
        TEST(indexer.local_index_to_global(180) == 540);
      }
    }

  } // end test with no order reduction

  //--------------------------------------------------------------------------//
  // TEST INDEXER WITH SPACE-ANGLE ORDER REDUCTION
  //--------------------------------------------------------------------------//

  {

    db->put<int>("erme_order_reduction", 3);
    ResponseIndexer indexer(db, nodes);
    TEST(indexer.number_nodes()                     == 4);
    TEST(indexer.number_node_moments(0)             == 108);
    TEST(indexer.response_index(0, 0, 0).azimuth    == 0);
    TEST(indexer.response_index(0, 0, 26).polar     == 2);
    TEST(indexer.response_index(0, 0, 26).azimuth   == 2);
    TEST(indexer.response_index(0, 0, 26).space0    == 0);

    // Loop over GLOBAL node indices
    for (int n = nodes->lower_bound(); n < nodes->upper_bound(); n++)
    {
      // Moment index on node
      int n_index = 0;
      // Surfaces
      for (int s = 0; s < 4; s++)
      {
        // Neighbor GLOBAL index and its surface
        int neigh_n = nodes->neighbor(n, s).neighbor();
        int neigh_s = nodes->neighbor(n, s).surface();

        // All moments on surface s
        for (int m = 0; m < indexer.number_surface_moments(n, s); m++, n_index++)
        {
          // If this surface is shared by a node, ensure we match orders
          if (nodes->neighbor(n, s).neighbor() >= 0)
          {
            ResponseIndex neigh_index = indexer.response_index(neigh_n, neigh_s, m);
            ResponseIndex index       = indexer.response_index(n, s, m);
            TEST(neigh_index.surface  == neigh_s);
            TEST(neigh_index.energy   == index.energy);
            TEST(neigh_index.polar    == index.polar);
            TEST(neigh_index.azimuth  == index.azimuth);
            TEST(neigh_index.space0   == index.space0);
            TEST(neigh_index.space1   == index.space1);
          }
          // Otherwise, ensure it's vacuum
          else
          {
            TEST(nodes->neighbor(n, s).neighbor() == erme_geometry::Node::VACUUM);
          }
        } // end moments
      } // surfaces
    } // nodes

  } // end test indexer with space-angle order reduction

  Comm::finalize();
  return 0;
}

/*
 *  This uses the smallest index set possible in order to help
 *  debug multiprocess indexing.
 */
int test_ResponseIndexer_zeroth_order(int argc, char *argv[])
{

  //--------------------------------------------------------------------------//
  // SETUP COMM
  //--------------------------------------------------------------------------//

  typedef serment_comm::Comm Comm;

  // Use one local group for 1 or 2 proc's, and 2 for more.
  Comm::initialize(argc, argv);
  if (Comm::size() < 3)
    Comm::setup_communicators(1);
  else
    Comm::setup_communicators(2);

  //--------------------------------------------------------------------------//
  // CREATE NODE LIST AND PARTITION
  //--------------------------------------------------------------------------//

  // Get the node list
  erme_geometry::NodeList::SP_nodelist nodes;
  if (Comm::rank() == 0)
    nodes = erme_geometry::cartesian_node_dummy_list_2d(0, 0, 0);

  // Partition the nodes
  erme_geometry::NodePartitioner partitioner;
  partitioner.partition(nodes);

  //--------------------------------------------------------------------------//
  // CREATE INDEXER
  //--------------------------------------------------------------------------//

  // Create parameter database
  ResponseIndexer::SP_db db(new detran_utilities::InputDB());
  db->put<int>("dimension", 2);

  db->put<int>("erme_order_reduction", 0);
  ResponseIndexer indexer(db, nodes);

  //--------------------------------------------------------------------------//
  // TEST
  //--------------------------------------------------------------------------//

  TEST(indexer.number_nodes()           == 4);
  TEST(indexer.number_node_moments(0)   == 4);
  TEST(indexer.number_global_moments()  == 16);

  if (Comm::size() < 3)
  {
    TEST(indexer.number_local_moments() == 16);
  }
  else
  {
    TEST(indexer.number_local_moments() == 8);
  }

  Comm::finalize();
  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_StateERME.cc
//----------------------------------------------------------------------------//
