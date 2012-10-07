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
   *   3.  one group diffusion slab evaluated at ten keffs
   *   4.  two group diffusion with second order space at three keffs
   *
   */

  // Load database
  ResponseDatabase::SP_rfdb db;
  db = ResponseDatabase::Create("test.h5");

  // Create response and fill
  SP_response rf(new NodeResponse(2, 2));
  db->get("node", rf, 0, 0.9); // left side
  db->get("node", rf, 1, 0.9); // right side

  // Test exact values.
  TEST(soft_equiv(rf->boundary_response(0, 0), 1.0));
  TEST(soft_equiv(rf->boundary_response(1, 0), 1.0));
  TEST(soft_equiv(rf->boundary_response(0, 1), 1.0));
  TEST(soft_equiv(rf->boundary_response(1, 1), 1.0));

  // Test interpolated values

  // Finish up
  Comm::finalize();
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_StateERME.cc
//---------------------------------------------------------------------------//
