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
#include "erme_geometry/NodeList.hh"
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

// Test of basic public interface
int test_NodeList(int argc, char *argv[])
{
  NodeList nodes = cartesian_node_detran_list_2d();
  TEST(nodes.number_global_nodes() == 4);
  TEST(nodes.neighbor(0, CartesianNode::TOP).neighbor()   == 2);
  TEST(nodes.neighbor(0, CartesianNode::TOP).surface()    == CartesianNode::BOTTOM);
  TEST(nodes.neighbor(1, CartesianNode::RIGHT).neighbor() == Node::VACUUM);
  TEST(nodes.node(2)->dimension() == 2);
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_NodeList.cc
//---------------------------------------------------------------------------//
