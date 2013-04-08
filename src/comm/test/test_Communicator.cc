//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_Communicator.cc
 *  @author Jeremy Roberts
 *  @date   Aug 19, 2012
 *  @brief  Test of communicator creation and use.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                       \
        FUNC(test_Communicator)

// Detran test
#include "utilities/TestDriver.hh"

#include "Comm.hh"

#include <cstdio>
#include <cmath>
#include <iostream>
#include <vector>

// Setup

using namespace serment_comm;
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

// Test of send and receive
int test_Communicator(int argc, char *argv[])
{

  // Initialize Comm
  Comm::initialize(argc, argv);

  // Set the number of local communicators
  unsigned int N = 0;
  if (Comm::size() == 1) N = 1;
  if (Comm::size() == 2) N = 2;
  if (Comm::size() == 3) N = 2;
  if (Comm::size() == 4) N = 2;
  if (Comm::size() == 5) N = 3;
  if (Comm::size() >= 6) return 0;

  int W = Comm::rank();

  // Initialize global and local
  Comm::setup_communicators(N);

  // Switch to local.
  Comm::set(local);
  int L = Comm::rank();

  // Switch to global, if applicable
  int G = -1;
  if (Comm::is_global())
  {
    Comm::set(global);
    G = Comm::rank();
  }

  // Switch to world.
  Comm::set(world);
  cout << "RANK = " << Comm::rank() << "W=" << W << " G=" << G << " L=" << L << endl;

  if (Comm::rank() == 0)
  {
    TEST(W == 0 and G == 0 and L == 0);
  }
  if (Comm::size() == 2)
  {
    if (Comm::rank() == 1)
    {
      TEST(W == 1 and G == 1 and L == 0);
    }
  }
  if (Comm::size() == 3)
  {
    if (Comm::rank() == 1)
    {
      TEST(W == 1 and G == -1 and L == 1);
    }
    if (Comm::rank() == 2)
    {
      TEST(W == 2 and G == 1 and L == 0);
    }
  }
  if (Comm::size() == 4)
  {
    if (Comm::rank() == 1)
    {
      TEST(W == 1 and G == -1 and L == 1);
    }
    if (Comm::rank() == 2)
    {
      TEST(W == 2 and G == 1 and L == 0);
    }
    if (Comm::rank() == 3)
    {
      TEST(W == 3 and G == -1 and L == 1);
    }
  }
  if (Comm::size() == 5)
  {
    if (Comm::rank() == 1)
    {
      TEST(W == 1 and G == -1 and L == 1);
    }
    if (Comm::rank() == 2)
    {
      TEST(W == 2 and G == 1 and L == 0);
    }
    if (Comm::rank() == 3)
    {
      TEST(W == 3 and G == -1 and L == 1);
    }
    if (Comm::rank() == 4)
    {
      TEST(W == 4 and G == 2 and L == 0);
    }
  }

  // Testing broadcast
  int data = 0;
  if (Comm::rank() == 0)
  {
    data = 1;
  }
  Comm::broadcast(&data, 1, 0);
  TEST(data == 1);

  cout << " RANK " << Comm::rank() << " is done." << endl;

  // Free the comm.
  Comm::free();

  // Finalize Comm
  Comm::finalize();

  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_Communicator.cc
//---------------------------------------------------------------------------//
