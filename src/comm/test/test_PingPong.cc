//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_Comm.cc
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 * \brief  Test of Vector class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                       \
        FUNC(test_PingPong_sendreceive) \
        FUNC(test_PingPong_bandwidth)   \
        FUNC(test_PingPong_latency)

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
int test_PingPong_sendreceive(int argc, char *argv[])
{

  // Initialize Comm
  Comm::initialize(argc, argv);

  // Test is valid for two processes only
  if (Comm::size() != 2) return 0;

  char   c = 0;
  int    i = 0;
  long   l = 0;
  float  f = 0;
  double d = 0;

  // assign on node 0
  if (Comm::rank() == 0)
  {
    c = 'A';
    i = 1;
    l = 1000;
    f = 1.5;
    d = 2.5;

    // send out data
    Comm::send(&c, 1, 1);
    Comm::send(&i, 1, 1);
    Comm::send(&l, 1, 1);
    Comm::send(&f, 1, 1);
    Comm::send(&d, 1, 1);

    // receive back
    Comm::receive(&c, 1, 1);
    Comm::receive(&i, 1, 1);
    Comm::receive(&l, 1, 1);
    Comm::receive(&f, 1, 1);
    Comm::receive(&d, 1, 1);

    // check values
    TEST(c == 'B');
    TEST(i == 2);
    TEST(l == 2000);
    TEST(soft_equiv(f, 2.5f));
    TEST(soft_equiv(d, 3.5));
  }

  // receive and send on node 1
  if (Comm::rank() == 1)
  {
    // receive from node 0
    Comm::receive(&c, 1, 0);
    Comm::receive(&i, 1, 0);
    Comm::receive(&l, 1, 0);
    Comm::receive(&f, 1, 0);
    Comm::receive(&d, 1, 0);

    // check values
    TEST(c == 'A');
    TEST(i == 1);
    TEST(l == 1000);
    TEST(soft_equiv(f, 1.5f));
    TEST(soft_equiv(d, 2.5));

    // assign new values
    c = 'B';
    i = 2;
    l = 2000;
    f = 2.5;
    d = 3.5;

    // send them back
    Comm::send(&c, 1, 0);
    Comm::send(&i, 1, 0);
    Comm::send(&l, 1, 0);
    Comm::send(&f, 1, 0);
    Comm::send(&d, 1, 0);

  }

  // Finalize Comm
  Comm::finalize();
  return 0;
}

// This tests the latency of the system.  The user must
// know how to setup the run to get processes on the
// nodes of interest.
int test_PingPong_bandwidth(int argc, char *argv[])
{

  // Initialize Comm
  Comm::initialize(argc, argv);

  // Test is valid for two processes only
  if (Comm::size() != 2) return 0;

  // Do the test for several array sizes of powers of 2
  for (int n = -1; n < 20; n++)
  {

    int count = 0;
    if (n >= 0) count = std::pow(2.0, n);

    // Create data
    std::vector<double> buffer(count, 1.0);

    // Start timer
    Comm::tic();

    // Do several tests.
    int number_trials = 100;
    for (int trial = 0; trial < number_trials; trial++)
    {

      // assign on node 0
      if (Comm::rank() == 0)
      {
        // send out data
        Comm::send(&buffer[0], count, 1);

        // receive back
        Comm::receive(&buffer[0], count, 1);
      }

      // receive and send on node 1
      if (Comm::rank() == 1)
      {

        // create data
        std::vector<double> buffer(count, 1.0);

        // receive from node 0
        Comm::receive(&buffer[0], count, 0);

        // send data
        Comm::send(&buffer[0], count, 0);
      }

    } // end trials

    Comm::global_barrier();

    // Print in a nice format
    if (Comm::rank() == 0)
    {
      double etime = Comm::toc();
      // 2 way * 64 bit/double * 1.0e6 bit/MB
      double ms = (double) count * 128 / 1.0e6;
      double bw = (double) number_trials * ms / etime;
      cout << ms << "  " << bw << endl;
    }

  } // end sizes

  // Finalize Comm
  Comm::finalize();
  return 0;
}


// This tests the latency of the system.  The user must
// know how to setup the run to get processes on the
// nodes of interest.
int test_PingPong_latency(int argc, char *argv[])
{

  // Initialize Comm
  Comm::initialize(argc, argv);

  // Test is valid for two processes only
  if (Comm::size() != 2) return 0;
  Comm::global_barrier();

  unsigned char msg = 'x';
  int reps = 500;

  if (Comm::rank() == 0)
  {
     // round-trip latency timing test
     std::printf("task %d has started...\n", Comm::rank());
     std::printf("Beginning latency timing test. Number of reps = %d.\n", reps);
     std::printf("***************************************************\n");
     std::printf("Rep#       T1               T2            deltaT\n");

     double time = 0.0;
     int average_time = 0.0;
     double total_time = 0.0;

     for (int n = 1; n <= reps; n++)
     {

       // Start timer
       Comm::tic();

       // Send message to worker
       Comm::send(&msg, 1, 1);

       // Wait to receive the echo reply from the worker
       Comm::receive(&msg, 1, 1);

       // Get elapsed time and add to total
       time = Comm::toc();
       if (n % 10 == 0) std::printf("%4d  %2.8f\n", n, time);
       total_time += time;

     }

     average_time = (total_time * 1.0e6) / reps;
     std::printf("***************************************************\n");
     std::printf("\n*** Avg round trip time = %d microseconds\n", average_time);
     std::printf("*** Avg one way latency = %d microseconds\n", average_time/2);

  }

  else if (Comm::rank() == 1)
  {
    std::printf("task %d has started...\n", Comm::rank());
    for (int n = 1; n <= reps; n++)
    {
       Comm::receive(&msg, 1, 0);
       Comm::send(&msg, 1, 0);
    }
  }

  // Finalize Comm
  Comm::finalize();
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_PingPong.cc
//---------------------------------------------------------------------------//
