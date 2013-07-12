//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_GlobalSolverPicard.cc
 *  @brief Test of GlobalSolverPicard class
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_GlobalSolverPicard)

#include "utilities/TestDriver.hh"
#include "erme_solver/ManagerERME.hh"
#include "erme_geometry/test/nodelist_fixture.hh"
#include <iostream>

using namespace serment_comm;
using namespace erme_solver;
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

int test_GlobalSolverPicard(int argc, char *argv[])
{
  // Parameter database
  ManagerERME::SP_db db = detran_utilities::InputDB::Create();

  // Setup the manager and comm
  ManagerERME manager(argc, argv);
  int ng = Comm::size() == 1 ? 1 : 2;
  db->put<std::string>("erme_solver_type", "picard");
  db->put<int>("comm_local_groups", 1);
  db->put<int>("dimension", 1);
  db->put<int>("erme_maximum_iterations", 10);
  db->put("basis_p_type", "jacobi");
  manager.build_comm(db);

  // Get nodes, build problem, and solve
  NodeList::SP_nodelist nodes = cartesian_node_detran_list_1d(0);//, 0, 0);
  manager.build_erme(nodes);
//  nodes->display();
  manager.solve();

  serment_comm::Comm::global_barrier();
  return 0;
}

//----------------------------------------------------------------------------//
//              end of file test_GlobalSolverPicard.cc
//----------------------------------------------------------------------------//
