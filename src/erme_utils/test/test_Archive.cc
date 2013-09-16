//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Archive.cc
 *  @brief test_Archive
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST           \
        FUNC(test_Archive)

#include "utilities/TestDriver.hh"
#include "Archive.hh"
#include "erme_geometry/test/nodelist_fixture.hh"

using namespace erme_utils;
using namespace erme_geometry;
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

int test_Archive(int argc, char *argv[])
{
  std::string filename = "test.archive";
  int number_nodes = 0;

  // Write
  if(0){
    Archive::SP_db db;
    Archive::SP_nodelist nodes;
    db = new detran_utilities::InputDB("testdb");
    db->put<int>("ip", 123);
    db->put<double>("dp", 0.456);
    //std::cout << db << std::endl;
    nodes = cartesian_node_dummy_list_1d();
    number_nodes = nodes->number_global_nodes();
    Archive A;
    A.archive(db, nodes, filename);
  }

  // Read
  if(1){
    Archive::SP_db db;
    Archive::SP_nodelist nodes;
    Archive A;
    A.unarchive(filename);
    db = A.get_db();
    nodes = A.get_nodes();
    TEST(db->get<int>("ip") == 123);
    TEST(soft_equiv(db->get<double>("dp"), 0.456));
    TEST(number_nodes = nodes->number_global_nodes());
  }

  return 0;
}
