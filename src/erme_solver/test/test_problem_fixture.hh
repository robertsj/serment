//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_problem_fixture.hh
 *  @brief Definition of several Detran example problems for testing
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef TEST_PROBLEM_FIXTURE_HH_
#define TEST_PROBLEM_FIXTURE_HH_

#include "erme_geometry/test/nodelist_fixture.hh"


/// Produce a test detran problem
static void erme_test_problem()
{

  // Get the detran nodelist fixture
  erme_geometry::NodeList::SP_nodelist nodes =
    erme_geometry::cartesian_node_detran_list_1d();

  //

}




#endif /* TEST_PROBLEM_FIXTURE_HH_ */


//----------------------------------------------------------------------------//
//              end of file test_problem_fixture.hh
//----------------------------------------------------------------------------//
