//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_GlobalSolverNewton.cc
 *  @brief Test of GlobalSolverNewton class
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_GlobalSolverNewton)

#include "utilities/TestDriver.hh"
#include "erme_solver/ManagerERME.hh"
#include "erme_geometry/test/nodelist_fixture.hh"
#include <iostream>

using namespace serment_comm;
using namespace erme_solver;
using namespace erme_geometry;
using namespace erme_response;
using namespace detran_test;
using namespace detran_utilities;
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

int test_GlobalSolverNewton(int argc, char *argv[])
{
  // Parameter database
  ManagerERME::SP_db db = detran_utilities::InputDB::Create();

  // Setup the manager and comm
  ManagerERME manager(argc, argv);
  int ng = Comm::size() == 1 ? 1 : 2;
  db->put<std::string>("erme_solver_type", "newton");
  db->put<int>("comm_local_groups", 1);
  db->put<int>("dimension", 1);
  db->put<int>("erme_maximum_iterations", 10);
  db->put<double>("erme_tolerance", 1.0e-12);
  manager.build_comm(db);

  // Get nodes, build problem, and solve
  NodeList::SP_nodelist nodes = cartesian_node_detran_list_1d(1);
  manager.build_erme(nodes);
  manager.solve();

  TEST(soft_equiv(manager.get_keff(), 0.996181414, 1.0e-8));

  serment_comm::Comm::global_barrier();
  return 0;
}

//----------------------------------------------------------------------------//
//              end of file test_GlobalSolverNewton.cc
//----------------------------------------------------------------------------//
