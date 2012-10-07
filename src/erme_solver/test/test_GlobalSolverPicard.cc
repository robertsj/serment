//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_GlobalSolverPicard.cc
 *  @brief  test_GlobalSolverPicard
 *  @author Jeremy Roberts
 *  @date   Oct 6, 2012
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_GlobalSolverPicard)

#include "utilities/TestDriver.hh"
#include "erme_response/NodeResponse.hh"
#include <iostream>

// Setup

using namespace erme_response;
using namespace detran_test;
using detran_utilities::soft_equiv;
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
 *  This test is homogeneous slab with
 *  vacuum boundaries based on database
 *  response functions.
 */

int test_GlobalSolverPicard(int argc, char *argv[])
{


}

//---------------------------------------------------------------------------//
//              end of file test_GlobalSolverPicard.cc
//---------------------------------------------------------------------------//
