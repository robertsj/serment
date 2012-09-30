//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_GlobalReduction.cc
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 * \brief  Test of Vector class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                  \
        FUNC(test_GlobalReduction)

// Detran test
#include "utilities/TestDriver.hh"

#include "Comm.hh"

#include <cstdio>
#include <cmath>
#include <iostream>
#include <vector>
#include <mpi.h>
// Setup

using namespace serment_comm;
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

// Test of send and receive
int test_GlobalReduction(int argc, char *argv[])
{

  // Initialize Comm
  Comm::initialize(argc, argv);

  //----------------------------------------------//
  // SCALAR REDUCTIONS
  //----------------------------------------------//

  // sum
  double value = Comm::rank();
  double ref = 0.5 * (Comm::size() * (Comm::size() - 1));
  Comm::global_sum(value);
  TEST(soft_equiv(value, ref));

  // product
  value = Comm::rank() + 1.0;
  ref = 1;
  for (int i = 1; i <= Comm::size(); i++)
    ref *= (double) i;
  Comm::global_prod(value);
  TEST(soft_equiv(value, ref));

  // min
  value = Comm::rank() + 1.0;
  ref = 1.0;
  Comm::global_min(value);
  TEST(soft_equiv(value, ref));

  // max
  value = Comm::rank() + 1.0;
  ref = (double) Comm::size();
  Comm::global_max(value);
  TEST(soft_equiv(value, ref));

  //----------------------------------------------//
  // ARRAY REDUCTIONS
  //----------------------------------------------//

  // sum
  double values[] = {Comm::rank(), Comm::rank(), Comm::rank()};
  ref = 0.5 * (Comm::size() * (Comm::size() - 1));
  Comm::global_sum(values, 3);
  TEST(soft_equiv(values[0], ref));
  TEST(soft_equiv(values[1], ref));
  TEST(soft_equiv(values[2], ref));

  // product
  values[0] = Comm::rank() + 1.0;
  values[1] = Comm::rank() + 1.0;
  values[2] = Comm::rank() + 1.0;
  ref = 1;
  for (int i = 1; i <= Comm::size(); i++)
    ref *= (double) i;
  Comm::global_prod(values, 3);
  TEST(soft_equiv(values[0], ref));
  TEST(soft_equiv(values[1], ref));
  TEST(soft_equiv(values[2], ref));

  // min
  values[0] = Comm::rank() + 1.0;
  values[1] = Comm::rank() + 1.0;
  values[2] = Comm::rank() + 1.0;
  ref = 1.0;
  Comm::global_min(values, 3);
  TEST(soft_equiv(values[0], ref));
  TEST(soft_equiv(values[1], ref));
  TEST(soft_equiv(values[2], ref));

  // max
  values[0] = Comm::rank() + 1.0;
  values[1] = Comm::rank() + 1.0;
  values[2] = Comm::rank() + 1.0;
  ref = (double) Comm::size();
  Comm::global_max(values, 3);
  TEST(soft_equiv(values[0], ref));
  TEST(soft_equiv(values[1], ref));
  TEST(soft_equiv(values[2], ref));

  // Finalize Comm
  Comm::finalize();
  return 0;
}


//---------------------------------------------------------------------------//
//              end of test_PingPong.cc
//---------------------------------------------------------------------------//
