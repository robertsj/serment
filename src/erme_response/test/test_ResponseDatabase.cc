//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_ResponseDatabase.cc
 *  @brief Test of ResponseDatabase class.
 *  @note  Copyright (C) 2013 Jeremy Roberts
 *  @todo  Add switch between R~k, R~1/k, and 1/R~1/k
 */
//----------------------------------------------------------------------------//

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

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

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
   * The database has four nodes from the 2-D IAEA benchmark.  The fourth
   * node is moderator and is evaluated once.  The others are fuel and
   * are evaluated for 4 keff values.  The responses have two groups and
   * are first order in space.
   */

  // Load database
  ResponseDatabase::SP_rfdb db;
  db = ResponseDatabase::Create("test.h5");

  SP_response rf(new NodeResponse(16, 4, 0));
  db->get("node3", rf, ResponseIndex(0,0,0,0,0,0,0,0,0), 1.0);
  TEST(soft_equiv(rf->boundary_response(0, 0), 0.137712, 1e-5));
  db->get("node3", rf, ResponseIndex(0,0,0,0,0,0,0,0,0), 1.3);
  TEST(soft_equiv(rf->boundary_response(0, 0), 0.137712, 1e-5));


  // Finish up
  Comm::finalize();
  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_StateERME.cc
//----------------------------------------------------------------------------//
