//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Jacobian.cc
 *  @brief Test of Jacobian and FullJacobian classes
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_Jacobian)           \
        FUNC(test_FullJacobian)

#include "utilities/TestDriver.hh"
#include "erme_solver/ManagerERME.hh"
#include "erme_solver/Jacobian.hh"
#include "erme_solver/FullJacobian.hh"
#include "erme_geometry/test/nodelist_fixture.hh"
#include <iostream>

using namespace serment_comm;
using namespace linear_algebra;
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

void update(ManagerERME &manager, double keff, std::string name)
{
  manager.get_server()->update(keff);
  Comm::set(global);
  if (Comm::is_global())
  {
    const int BINARY = manager.get_responses()->R->BINARY;
    manager.get_responses()->R->update();
    manager.get_responses()->F->update();
    manager.get_responses()->A->update();
    manager.get_responses()->L->update();
    manager.get_responses()->R->display(BINARY, "jac_R"+name+".out");
    manager.get_responses()->M->display(BINARY, "jac_M"+name+".out");
    manager.get_responses()->F->display(BINARY, "jac_F"+name+".out");
    manager.get_responses()->A->display(BINARY, "jac_A"+name+".out");
    manager.get_responses()->L->display(BINARY, "jac_L"+name+".out");
    manager.get_responses()->L->leakage_vector().display(BINARY, "jac_Lv"+name+".out");
  }
  Comm::set(world);
}

int test_Jacobian(int argc, char *argv[])
{
  // Parameter database
  ManagerERME::SP_db db = detran_utilities::InputDB::Create();

  // Setup the manager and comm
  ManagerERME manager(argc, argv);
  int ng = Comm::size() == 1 ? 1 : 2;
  db->put<int>("comm_local_groups", ng);
  db->put<int>("dimension", 1);
  manager.build_comm(db);

  // Get nodes, build problem, and retrieve responses
  NodeList::SP_nodelist nodes = cartesian_node_dummy_list_1d();
  manager.build_erme(nodes);

  ManagerERME::SP_responsecontainer r = manager.get_responses();

  update(manager, 1.0001, "1");
  update(manager, 1.0000, "0");

  Jacobian jacobian(manager.get_server(), manager.get_indexer(), r, 0.0001);

  // Jacobian size
  int m = 0;
  // Unknowns
  Vector::SP_vector x;

  if (Comm::is_global())
  {
    Comm::set(global);
    cout << " one " << endl;
    m = r->R->number_local_rows();
    if (Comm::is_last()) m += 2;
    x = new Vector(m, 1.0);
    Comm::set(world);
  }

  cout << " two " << endl;

  jacobian.update(x);
  jacobian.matrix()->display(r->L->BINARY, "jac.out");


  cout << " three " << endl;

  linear_algebra::Vector jacobian_times_x(m, 0.0);
  jacobian.multiply(*x, jacobian_times_x);
  jacobian_times_x.display(r->A->BINARY, "jac_x.out");

  double ref[] = {-2.000000000000000e+00, -2.000000000000000e+00,
      -2.000000000000000e+00, 3.599980000600156e+05, 3.599980000600156e+05,
      3.599980000600156e+05, 3.599980000600156e+05, 3.599980000600156e+05,
      3.599980000600156e+05, 3.599980000600156e+05, 3.599980000600156e+05,
      3.599980000600156e+05, 3.599980000600156e+05, 3.599980000600156e+05,
      3.599980000600156e+05, 3.599980000600156e+05, 3.599980000600156e+05,
      3.599980000600156e+05, 3.599980000600156e+05, 3.599980000600156e+05,
      3.599980000600156e+05, -2.000000000000000e+00, -2.000000000000000e+00,
      -2.000000000000000e+00, -2.880024119519948e+06, -2.400000000000000e+01};

  for (int i = 0; i < jacobian_times_x.local_size(); ++i)
  {
    TEST(soft_equiv(jacobian_times_x[i], ref[i]));
  }
  Comm::global_barrier();
  return 0;
}

int test_FullJacobian(int argc, char *argv[])
{
  // Parameter database
  ManagerERME::SP_db db = detran_utilities::InputDB::Create();

  // Setup the manager and comm
  ManagerERME manager(argc, argv);
  int ng = Comm::size() == 1 ? 1 : 2;
  db->put<int>("comm_local_groups", ng);
  db->put<int>("dimension", 1);
  manager.build_comm(db);

  // Get nodes, build problem, and retrieve responses
  NodeList::SP_nodelist nodes = cartesian_node_dummy_list_1d();
  manager.build_erme(nodes);

  ManagerERME::SP_responsecontainer r = manager.get_responses();

  update(manager, 1.0001, "1");
  update(manager, 1.0000, "0");

  FullJacobian jacobian(manager.get_server(), manager.get_indexer(), r, 0.0001);
//  std::cout << jacobian.matrix()->number_global_rows() << std::endl;
//  return 0;
  // Jacobian size
  int m = 0;
  // Unknowns
  Vector::SP_vector x;

  if (Comm::is_global())
  {
    Comm::set(global);
    cout << " one " << endl;
    m = r->R->number_local_rows();
    if (Comm::is_last()) m += 2;
    x = new Vector(m, 1.0);
    Comm::set(world);
  }

  cout << " two " << endl;

  jacobian.update(x);
  jacobian.matrix()->display(r->L->BINARY, "fulljac.out");

  cout << " three " << endl;

  linear_algebra::Vector jacobian_times_x(m, 0.0);
  x->display(x->BINARY, "v_in_a.out");
  jacobian.multiply(*x, jacobian_times_x);
  jacobian_times_x.display(r->A->BINARY, "fulljac_x.out");

  double ref[] = {-2.000000000000000e+00, -2.000000000000000e+00,
      -2.000000000000000e+00, 3.599980000600156e+05, 3.599980000600156e+05,
      3.599980000600156e+05, 3.599980000600156e+05, 3.599980000600156e+05,
      3.599980000600156e+05, 3.599980000600156e+05, 3.599980000600156e+05,
      3.599980000600156e+05, 3.599980000600156e+05, 3.599980000600156e+05,
      3.599980000600156e+05, 3.599980000600156e+05, 3.599980000600156e+05,
      3.599980000600156e+05, 3.599980000600156e+05, 3.599980000600156e+05,
      3.599980000600156e+05, -2.000000000000000e+00, -2.000000000000000e+00,
      -2.000000000000000e+00, -2.880024119519948e+06, -2.400000000000000e+01};

  for (int i = 0; i < 5; ++i)
  {
    TEST(soft_equiv(jacobian_times_x[i], ref[i]));
  }

  Comm::global_barrier();
  return 0;
}

//----------------------------------------------------------------------------//
//              end of file test_Jacobian.cc
//----------------------------------------------------------------------------//
