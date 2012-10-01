//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_ResponseDatabase.cc
 * \author Jeremy Roberts
 * \date   Oct 1, 2012
 * \brief  Test of ResponseDatabase class.
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
 *
 */

int test_ResponseDatabase(int argc, char *argv[])
{

  typedef NodeResponse::SP_response SP_response;

  //-------------------------------------------------------------------------//
  // SETUP COMM
  //-------------------------------------------------------------------------//

  typedef serment_comm::Comm Comm;

  // Initialize comm
  Comm::initialize(argc, argv);

  // Setup local and global communicators.
  int number_local_comm = 1;
  if (Comm::size() > 2) number_local_comm = 2;
  Comm::setup_communicators(number_local_comm);

//  // Create server
  ResponseDatabase::SP_rfdb db;
  db = ResponseDatabase::Create("test.h5");
//  cout << db->filename() << endl;

  // Finish up
  Comm::finalize();
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_StateERME.cc
//---------------------------------------------------------------------------//
