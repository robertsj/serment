//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_ResponseDatabase.cc
 *  @author Jeremy Roberts
 *  @date   Oct 1, 2012
 *  @brief  Test of ResponseDatabase class.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                   \
        FUNC(test_ResponseDatabase)

#include "utilities/TestDriver.hh"
#include "erme_response/ResponseDatabase.hh"
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
 *  This reads from a test HDF5 database file.
 */

int test_ResponseDatabase(int argc, char *argv[])
{

  typedef NodeResponse::SP_response SP_response;
  typedef serment_comm::Comm Comm;

  // Initialize comm
  Comm::initialize(argc, argv);

  // Setup local and global communicators.
  int number_local_comm = 1;
  if (Comm::size() > 2) number_local_comm = 2;
  Comm::setup_communicators(number_local_comm);

  /*
   * The database has four nodes:
   *
   *   1.  one group diffusion slab evaluated at one keff
   *   2.  one group diffusion slab evaluated at two keffs
   *   3.  one group diffusion slab evaluated at 5 keffs
   *   4.  two group diffusion with second order space at 5 keffs
   *
   */

  // Load database
  ResponseDatabase::SP_rfdb db;
  db = ResponseDatabase::Create("test.h5");

  // case 1
  {
    // Create response and fill
    SP_response rf(new NodeResponse(2, 2));
    db->get("node1", rf, ResponseIndex(0,0,0,0,0,0,0,0,0), 1.0); // left side
    db->get("node1", rf, ResponseIndex(0,0,0,0,0,0,0,0,1), 1.0); // right side
    // Test exact values.
    TEST(soft_equiv(rf->boundary_response(0, 0), 0.9375));
    TEST(soft_equiv(rf->boundary_response(1, 0), 0.0625));
    TEST(soft_equiv(rf->boundary_response(0, 1), 0.0625));
    TEST(soft_equiv(rf->boundary_response(1, 1), 0.9375));
  }

  // case 2
  {
    // Create response and fill
    SP_response rf(new NodeResponse(2, 2));
    db->get("node2", rf, ResponseIndex(0,0,0,0,0,0,0,0,0), 1.0); // left side
    db->get("node2", rf, ResponseIndex(0,0,0,0,0,0,0,0,1), 1.0); // right side

    // Test interpolated values
    double ref0 = 0.5 * (1.521990 + 0.610188);
    double ref1 = 0.5 * (0.842184 + 0.0004057);
    TEST(soft_equiv(rf->boundary_response(0, 0), ref0, 1e-5));
    TEST(soft_equiv(rf->boundary_response(1, 0), ref1, 1e-5));
    TEST(soft_equiv(rf->boundary_response(0, 1), ref1, 1e-5));
    TEST(soft_equiv(rf->boundary_response(1, 1), ref0, 1e-5));
  }

  // case 3
  {
    // Create response and fill
    SP_response rf(new NodeResponse(2, 2));
    db->get("node3", rf, ResponseIndex(0,0,0,0,0,0,0,0,0), 1.05); // left side
    db->get("node3", rf, ResponseIndex(0,0,0,0,0,0,0,0,1), 1.05); // right side

    // Test interpolated values
    double ref0 = 0.5 * (0.9375 + 0.6101880);
    double ref1 = 0.5 * (0.0625 + 0.0004057);
    TEST(soft_equiv(rf->boundary_response(0, 0), ref0, 1e-5));
    TEST(soft_equiv(rf->boundary_response(1, 0), ref1, 1e-5));
    TEST(soft_equiv(rf->boundary_response(0, 1), ref1, 1e-5));
    TEST(soft_equiv(rf->boundary_response(1, 1), ref0, 1e-5));
  }

  // case 4
  {
    // Create response and fill
    SP_response rf(new NodeResponse(16, 4));
    for (int i = 0; i < 16; ++i)
      db->get("node4", rf, ResponseIndex(0,0,0,0,0,0,0,0,i), 1.05);

    // Test interpolated values
    double ref0 = 0.5 * (0.362366  + 0.349561);
    double ref1 = 0;
    double ref2 = 6.18677641e-02;
    TEST(soft_equiv(rf->boundary_response(0, 0), ref0, 1e-5));
    TEST(soft_equiv(rf->boundary_response(1, 0), ref1, 1e-5));
    TEST(soft_equiv(rf->boundary_response(0, 1), ref1, 1e-5));
    TEST(soft_equiv(rf->boundary_response(1, 1), ref2, 1e-5));
  }

  // Finish up
  Comm::finalize();
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_StateERME.cc
//---------------------------------------------------------------------------//
