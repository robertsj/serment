//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_Interpolation.cc
 *  @author Jeremy Roberts
 *  @date   Oct 1, 2012
 *  @brief  Test of Interpolation class.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                   \
        FUNC(test_Interpolation)

#include "utilities/TestDriver.hh"
#include "erme_response/Interpolation.hh"
#include <iostream>

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


int test_Interpolation(int argc, char *argv[])
{
  // f(x) = x^3
  double x0 = 0.0; double f0 = 0.0;
  double x1 = 1.0; double f1 = 1.0;
  double x2 = 1.5; double f2 = 3.375;
  double x3 = 2.0; double f3 = 8.0;

  // test linear
  double f_test = interpolate_linear(0.5, x0, x1, f0, f1);
  TEST(soft_equiv(f_test, 0.5));

  // test quadratic
  double f_test = interpolate_quadratic(0.5, x0, x1, x2, f0, f1, f2);
  TEST(soft_equiv(f_test, 0.5));

  // test linear
  double f_test = interpolate_cubic(0.5, x0, x1, x2, x3, f0, f1, f2, f3);
  TEST(soft_equiv(f_test, 0.5));

  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_StateERME.cc
//---------------------------------------------------------------------------//
