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

#include "TestDriver.hh"
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
  typedef serment_comm::Comm Comm;

  // Initialize Comm
  Comm::initialize(argc, argv);

  // Nodes
  NodeList nodes;
  if (Comm::rank() == 0) nodes = cartesian_node_detran_list_2d();

  // Create partitioner
  NodePartitioner partitioner;

  // Partition nodes
  partitioner.partition(nodes);

  // TEST
  if (Comm::size() == 1)
  {
    TEST(nodes.number_nodes() == 4);
  }
  else if (Comm::size() == 2)
  {
    TEST(nodes.number_nodes() == 2);
  }
  else if (Comm::size() == 3)
  {
    if (Comm::rank() == 2)
    {
      TEST(nodes.number_nodes() == 2);
    }
    else
    {
      TEST(nodes.number_nodes() == 1);
    }
  }
  else
  {
    TEST(nodes.number_nodes() == 1);
  }

  // Finalize Comm
  Comm::finalize();
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_NodePartitioner.cc
//---------------------------------------------------------------------------//
