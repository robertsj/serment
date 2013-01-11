//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_ResponseSourceDetran.cc
 *  @author Jeremy Roberts
 *  @date   Aug 19, 2012
 *  @brief  Test of ResponseServer class.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST               \
        FUNC(test_ResponseSourceDetran)

#include "utilities/TestDriver.hh"
#include "erme_response/ResponseIndexer.hh"
#include "erme_response/ResponseServer.hh"
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

int test_ResponseSourceDetran(int argc, char *argv[])
{

  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_StateERME.cc
//---------------------------------------------------------------------------//
