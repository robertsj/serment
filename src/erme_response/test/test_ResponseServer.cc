//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_ResponseServer.cc
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 * \brief  Test of ResponseServer class.
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

/*
 *  This tests how the response server is used within
 *  response operators.  In general, the program flow is as
 *  follows:
 *
 *    world   setup problem
 *            partition nodes
 *            setup indexer
 *            setup problem
 *            broadcast initial keff guess
 *            server->update(initial_guess)
 *    local   compute responses
 *    global  do petsc stuff
 *    world   broadcast keff and so on
 */

int test_ResponseServer(int argc, char *argv[])
{

  typedef NodeResponse::SP_response SP_response;

  using detran::soft_equiv;

  //-------------------------------------------------------------------------//
  // SETUP COMM
  //-------------------------------------------------------------------------//

  typedef serment_comm::Comm Comm;

  // Initialize comm
  Comm::initialize(argc, argv);

  // Setup local and global communicators.
  int number_local_comm = 1;
  if (Comm::size() > 2) number_local_comm = 2;
  Comm::setup_communicators(number_local_comm);

  //-------------------------------------------------------------------------//
  // WORLD
  //-------------------------------------------------------------------------//

  // Get the node list
  erme_geometry::NodeList::SP_nodelist nodes;
  if (Comm::rank() == 0)
    nodes = erme_geometry::cartesian_node_dummy_list_2d(0, 0, 0);

  // Partition the nodes
  erme_geometry::NodePartitioner partitioner;
  partitioner.partition(nodes);

  // Create parameter database
  ResponseIndexer::SP_db db(new detran::InputDB());
  db->put<int>("dimension", 2);

  // Create indexer
  db->put<int>("erme_order_reduction", 3);
  ResponseIndexer::SP_indexer indexer(new ResponseIndexer(db, nodes));

  // Create server
  ResponseServer server(nodes, indexer);

  // Update
  server.update(1.0);
  Comm::global_barrier();

  //-------------------------------------------------------------------------//
  // LOCAL
  //-------------------------------------------------------------------------//

  // Switch to local
  Comm::set(serment_comm::local);

  // Only root process has data at this point
  if (Comm::world_rank() == 0)
  {

  for (int n = nodes->lower_bound(); n < nodes->upper_bound(); n++)
  {
    // Response for this node
    SP_response r = server.response(nodes->local_index(n));

    // Loop over surfaces
    for (int s = 0; s < nodes->node(n)->number_surfaces(); s++)
    {
      // Loop over surface moments
      for (int m = 0; m < indexer->number_surface_moments(n, s); m++)
      {
        ResponseIndex index = indexer->response_index(n, s, m);
        unsigned int in = index.local;
        double value = 1000000.0 * index.node +
                        100000.0 * index.surface +
                         10000.0 * index.polar +
                          1000.0 * index.azimuth +
                           100.0 * index.space0 +
                            10.0 * index.space1 +
                             1.0 * index.energy;
        for (int out = 0; out < r->size(); out++)
        {
          //cout << " out = " << out << " br =  " << r->boundary_response(out, in) << " ref = " << value + 0.1 << endl;
          TEST(soft_equiv(r->boundary_response(out, in), value + 0.1));
        }
        TEST(soft_equiv(r->fission_response(in),    value + 0.2));
        TEST(soft_equiv(r->absorption_response(in), value + 0.3));
        for (int s = 0; s < r->number_surfaces(); s++)
        {
          TEST(soft_equiv(r->leakage_response(s, in), value + 0.4 + 0.01 * s));
        }
      }
    }

  } // end node loop

  } // end local root (aka global) block

  //-------------------------------------------------------------------------//
  // WORLD
  //-------------------------------------------------------------------------//

  Comm::set(serment_comm::world);
  Comm::finalize();
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_StateERME.cc
//---------------------------------------------------------------------------//
