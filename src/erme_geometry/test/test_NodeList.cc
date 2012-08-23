//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_NodeList.cc
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 * \brief  Test of NodeList class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST         \
        FUNC(test_NodeList)

#include "TestDriver.hh"
#include "NodeList.hh"
#include "node_fixture.hh"
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

// Test of basic public interface
int test_NodeList(int argc, char *argv[])
{

  // Create node list
  NodeList nodes;

  // Get 4 test nodes
  NodeList::SP_node node0 = serment_test::cartesian_node_detran(2);
  NodeList::SP_node node1 = serment_test::cartesian_node_detran(2);
  NodeList::SP_node node2 = serment_test::cartesian_node_detran(2);
  NodeList::SP_node node3 = serment_test::cartesian_node_detran(2);

  // Create neighbor lists.
  NodeList::vec_int neigh0(4, Node::VACUUM);
  NodeList::vec_int neigh1(4, Node::VACUUM);
  NodeList::vec_int neigh2(4, Node::VACUUM);
  NodeList::vec_int neigh3(4, Node::VACUUM);
  //  --- ---
  // | 2 | 3 |
  //  --- ---
  // | 0 | 1 |   with vacuum everywhere else
  //  --- ---
  neigh0[CartesianNode::TOP]    = 2;
  neigh0[CartesianNode::RIGHT]  = 1;
  neigh1[CartesianNode::LEFT]   = 0;
  neigh1[CartesianNode::TOP]    = 3;
  neigh2[CartesianNode::BOTTOM] = 0;
  neigh2[CartesianNode::RIGHT]  = 3;
  neigh3[CartesianNode::BOTTOM] = 1;
  neigh3[CartesianNode::LEFT]   = 2;

  // Add nodes
  nodes.add_node(node0, neigh0);
  nodes.add_node(node1, neigh1);
  nodes.add_node(node2, neigh2);
  nodes.add_node(node3, neigh3);
  nodes.finalize();

  TEST(nodes.number_nodes() == 4);
  TEST(nodes.neighbor(0, CartesianNode::TOP)   == 2);
  TEST(nodes.neighbor(1, CartesianNode::RIGHT) == Node::VACUUM);
  TEST(nodes.node(2)->dimension() == 2);
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_NodeList.cc
//---------------------------------------------------------------------------//
