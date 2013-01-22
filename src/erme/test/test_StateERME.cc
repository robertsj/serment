//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_StateERME.cc
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 * \brief  Test of StateERME class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST           \
        FUNC(test_StateERME)

#include "utilities/TestDriver.hh"
#include "StateERME.hh"
#include <iostream>

// Setup

using namespace erme;
using namespace detran_test;
using detran_utilities::soft_equiv;
using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
  RUN(argc, argv);
  PetscFinalize();
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

// Test of basic public interface
int test_StateERME(int argc, char *argv[])
{

  // Create node list
  StateERME state(10);

  linear_algebra::Vector v(10, 1.0);

  state.set_k(1.12);
  state.set_lambda(1.13);
  TEST(soft_equiv(state.k(),      1.12));
  TEST(soft_equiv(state.lambda(), 1.13));

  state.update(v, 2.0, 3.0);

  TEST(soft_equiv(state.k(),      2.0));
  TEST(soft_equiv(state.lambda(), 3.0));

  return 0;

}

//---------------------------------------------------------------------------//
//              end of test_StateERME.cc
//---------------------------------------------------------------------------//
