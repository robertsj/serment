//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_NodeResponse.cc
 *  @brief Test of NodeResponse class
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST           \
        FUNC(test_NodeResponse)

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

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

// Test of basic public interface
int test_NodeResponse(int argc, char *argv[])
{
  NodeResponse response(16, 4);
  typedef unsigned int size_t;

  for (size_t in = 0; in < 16; in++)
  {
    for (size_t out = 0; out < 16; out++)
      response.boundary_response(out, in) = in + out;
    response.fission_response(in)    = 1.0 * in;
    response.absorption_response(in) = 2.0 * in;
    for (size_t surface = 0; surface < 4; surface++)
      response.leakage_response(surface, in) = 100.0*surface + in;
  }

  double *b = &response.boundary_response(0, 5);
  double *f = &response.fission_response(0);
  double *a = &response.absorption_response(0);
  double *l = &response.leakage_response(0, 3);
  for (size_t i = 0; i < 16; i++)
  {
    TEST(soft_equiv(b[i], 5.0 + i));
    TEST(soft_equiv(f[i], 1.0 * i));
    TEST(soft_equiv(a[i], 2.0 * i));
  }
  TEST(soft_equiv(l[0],   3.0));
  TEST(soft_equiv(l[1], 103.0));
  TEST(soft_equiv(l[2], 203.0));
  TEST(soft_equiv(l[3], 303.0));

  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_StateERME.cc
//----------------------------------------------------------------------------//
