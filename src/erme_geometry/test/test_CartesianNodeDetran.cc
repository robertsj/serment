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
#define TEST_LIST                      \
        FUNC(test_CartesianNodeDetran)

#include "TestDriver.hh"
#include "CartesianNodeDetran.hh"
#include <iostream>
#include "node_fixture.hh"

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

  Node_T::SP_node node = serment_test::cartesian_node_detran(1);

  cout << "Color = " << node->color(Node_T::Point(0.5, 0, 0)) << endl;



  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_CartesianNodeDetran.cc
//---------------------------------------------------------------------------//
