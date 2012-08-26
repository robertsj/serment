//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_CartesianNodeDetran.cc
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 * \brief  Test of CartesianNodeDetran class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                               \
        FUNC(test_CartesianNodeDetran)          \
        FUNC(test_CartesianNodeDetran_serialize)

#include "TestDriver.hh"
#include "erme_geometry/CartesianNodeDetran.hh"
#include "erme_geometry/test/node_fixture.hh"

#include <iostream>
#include <fstream>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

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
int test_CartesianNodeDetran(int argc, char *argv[])
{
  typedef CartesianNodeDetran Node_T;

  Node_T::SP_node node = cartesian_node_detran(1);

  cout << "Color = " << node->color(Node_T::Point(0.5, 0, 0)) << endl;
  return 0;
}

// Test of basic public interface
int test_CartesianNodeDetran_serialize(int argc, char *argv[])
{

  typedef Node Node_T;

  // Get node and pack.
  {
    Node_T::SP_node node = cartesian_node_detran(1);
    std::ofstream ofs("node.archive");
    boost::archive::binary_oarchive oa(ofs);
    oa << node;
    ofs.close();
  }

  // Unpack
  {
    Node_T::SP_node node;
    std::ifstream ifs("node.archive");
    boost::archive::binary_iarchive ia(ifs);
    ia >> node;
    ifs.close();
  }

  return 0;
}


//---------------------------------------------------------------------------//
//              end of test_CartesianNodeDetran.cc
//---------------------------------------------------------------------------//
