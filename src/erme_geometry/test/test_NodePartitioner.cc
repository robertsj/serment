//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_NodePartitioner.cc
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 * \brief  Test of test_NodePartitioner class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                 \
        FUNC(test_NodePartitioner)

#include "utilities/TestDriver.hh"
#include "erme_geometry/NodePartitioner.hh"
#include "erme_geometry/test/nodelist_fixture.hh"
#include <iostream>

// Setup

using namespace erme_geometry;
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

int test_NodePartitioner(int argc, char *argv[])
{

  //-------------------------------------------------------------------------//
  // SETUP COMM
  //-------------------------------------------------------------------------//

  typedef serment_comm::Comm Comm;

  // Initialize Comm
  Comm::initialize(argc, argv);

  // Setup communicators.  We use one local communicator for 1 or 2
  // processes; otherwise, we use two local communicators.
  int number_local_comm = 1;
  if (Comm::size() > 2)
    number_local_comm = 2;
  Comm::setup_communicators(number_local_comm);

  // Ensure we've got the right local groups.
  if (Comm::size() < 3)
  {
    TEST(Comm::local_group() == 0);
  }
  else
  {
    // 3 procs  [0,1]   [2]
    // 4 procs  [0,1]   [2,3]
    // 5 procs  [0,1,2] [3,4]
    if (Comm::rank() < Comm::size()/2 + Comm::size() % 2)
    {
      TEST(Comm::local_group() == 0);
    }
    else
    {
      TEST(Comm::local_group() == 1);
    }
  }

  //-------------------------------------------------------------------------//
  // CREATE NODES AND PARTITION
  //-------------------------------------------------------------------------//

  // Nodes
  NodeList::SP_nodelist nodes;
  if (Comm::rank() == 0) nodes = cartesian_node_dummy_list_2d();

  // Create partitioner
  NodePartitioner partitioner;

  // Partition nodes
  partitioner.partition(nodes);
  // Must be back in world.
  TEST(serment_comm::communicator == serment_comm::world);

  //-------------------------------------------------------------------------//
  // TEST
  //-------------------------------------------------------------------------//

  if (Comm::size() < 3)
  {
    // We have one local communicator, so all nodes are local
    TEST(nodes->number_global_nodes() == 4);
    TEST(nodes->number_local_nodes()  == 4);
  }
  else
  {
    // We have two local communicators, so just half the nodes are local
    TEST(nodes->number_global_nodes() == 4);
    TEST(nodes->number_local_nodes()  == 2);
    if (Comm::local_group() == 0)
    {
      TEST(nodes->lower_bound() == 0);
      TEST(nodes->upper_bound() == 2);
    }
    else
    {
      TEST(nodes->lower_bound() == 2);
      TEST(nodes->upper_bound() == 4);
    }
  }

  // Finalize Comm
  Comm::finalize();
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_NodePartitioner.cc
//---------------------------------------------------------------------------//
