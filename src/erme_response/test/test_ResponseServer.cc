//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_ResponseServer.cc
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 * \brief  Test of ResponseIndexer class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST               \
        FUNC(test_ResponseServer)

#include "TestDriver.hh"
#include "erme_response/ResponseIndexer.hh"
#include "erme_response/ResponseServer.hh"
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
  typedef NodeResponse::SP_response SP_response;

  Comm::initialize(argc, argv);

  // Get the node list
  erme_geometry::NodeList nodes;
  if (Comm::rank() == 0)
    nodes = erme_geometry::cartesian_node_dummy_list_2d();

  // Partition the nodes
  erme_geometry::NodePartitioner partitioner;
  partitioner.partition(nodes);

  // Create parameter database
  ResponseIndexer::SP_db db(new detran::InputDB());
  db->put<int>("dimension", 2);

  // Create indexer
  db->put<int>("erme_order_reduction", 3);
  ResponseIndexer indexer(db, nodes);

  // Create server
  ResponseServer server(nodes, indexer);

  // Update
  server.update(1.0);

  // Switch to global and get responses for each node.
  Comm::set(global);
  for (int n = nodes.lower_bound(); n < nodes.upper_bound(); n++)
  {
    SP_response r = server.response(n);

    // Test the response.  Our access as used here is the same
    // as used to fill global matrices.
    for (int m = 0; m < indexer.number_node_moments(n); m++)
    {
      for (int mm = 0; mm < indexer.number_node_moments(n); mm++)
        TEST(detran::soft_equiv(r->boundary_response(m, mm), 1.0 * n));
      TEST(detran::soft_equiv(r->fission_response(m),    2.0 * n));
      TEST(detran::soft_equiv(r->absorption_response(m), 3.0 * n));
      TEST(detran::soft_equiv(r->leakage_response(0, m), 4.0 * n));
      TEST(detran::soft_equiv(r->leakage_response(1, m), 5.0 * n));
      TEST(detran::soft_equiv(r->leakage_response(2, m), 6.0 * n));
      TEST(detran::soft_equiv(r->leakage_response(3, m), 7.0 * n));
    }
  }
  Comm::set(world);

  Comm::finalize();
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_StateERME.cc
//---------------------------------------------------------------------------//
