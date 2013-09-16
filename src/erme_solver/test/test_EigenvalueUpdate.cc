//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_EigenvalueUpdate.cc
 *  @brief Test of GlobalSolverNewton class
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                        \
        FUNC(test_EigenvalueUpdate)      \
        FUNC(test_SteffensenUpdate)      \
        FUNC(test_AnghelUpdateEXP)       \
        FUNC(test_AnghelUpdateLINLIN)    \
        FUNC(test_AnghelUpdateINVLIN)    \
        FUNC(test_AnghelUpdateINVINV)


#include "utilities/TestDriver.hh"
#include "erme_solver/ManagerERME.hh"
#include "erme_solver/SteffensenUpdate.hh"
#include "erme_solver/AnghelUpdate.hh"
#include "erme_geometry/test/nodelist_fixture.hh"
#include <iostream>

using namespace serment_comm;
using namespace erme_solver;
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

ManagerERME::SP_manager manager;

//----------------------------------------------------------------------------//
void build_manager(int argc, char *argv[])
{
  ManagerERME::initialize(argc, argv);

  // Parameter database
  ManagerERME::SP_db db = detran_utilities::InputDB::Create();
  manager = new ManagerERME(argc, argv);
  int ng = Comm::size() == 1 ? 1 : 2;
  db->put<std::string>("erme_solver_type", "picard");
  db->put<int>("comm_local_groups", 1);
  manager->build_comm(db);
}

EigenvalueUpdate::SP_update update;
double keff[5] = {1.0000000000000000, 0.996399409597206, 0.996194383674861,
    0.996182187777422, 0.9961814143300799};
double lambda[4] = {0.9639030592301611, 0.9977908937167308, 0.9998680268529435,
    0.9999921276201872};

//----------------------------------------------------------------------------//
void build_update(std::string key, int anghel_scheme = 0)
{
  if (key == "default")
    update = new EigenvalueUpdate();
  else if (key == "steffensen")
    update = new SteffensenUpdate(3);
  else if (key == "anghel")
    update = new AnghelUpdate(anghel_scheme, false);
}

linear_algebra::Vector::SP_vector J;

//----------------------------------------------------------------------------//
int test_EigenvalueUpdate(int argc, char *argv[])
{
  {
    build_manager(argc, argv);
    build_update("default");
    double ref[4] = {1.0, 0.996399409597, 0.996194383675, 0.996182187777};
    for (int i = 0; i < 4; ++i)
    {
      TEST(soft_equiv(update->compute(keff[i], lambda[i], J), ref[i]));
    }
  }
  ManagerERME::finalize();
  return 0;
}

//----------------------------------------------------------------------------//
int test_SteffensenUpdate(int argc, char *argv[])
{
  {
    build_manager(argc, argv);
    build_update("steffensen");
    double ref[4] = {keff[0], keff[1], keff[2], 0.996181416424961};
    for (int i = 0; i < 4; ++i)
    {
      double val = update->compute(keff[i], lambda[i], J);
      printf("%16.12f %16.12f %16.12f \n", keff[i], val, ref[i]);
      TEST(soft_equiv(val, ref[i]));
    }
  }
  ManagerERME::finalize();
  return 0;
}

//----------------------------------------------------------------------------//
int test_AnghelUpdateEXP(int argc, char *argv[])
{
  {
    build_manager(argc, argv);
    build_update("anghel", AnghelUpdate::EXP);
    double ref[4] = {keff[0], keff[1], 0.996169838031401, 0.996181374348308};
    for (int i = 1; i < 4; ++i)
    {
      double val = update->compute(keff[i], lambda[i-1], J);
      printf("%16.12f %16.12f %16.12f \n", keff[i], val, ref[i]);
      TEST(soft_equiv(val, ref[i]));
    }
  }
  ManagerERME::finalize();
  return 0;
}

//----------------------------------------------------------------------------//
int test_AnghelUpdateLINLIN(int argc, char *argv[])
{
  {
    build_manager(argc, argv);
    build_update("anghel", AnghelUpdate::LINLIN);
    double ref[4] = {keff[0], keff[1], 0.996164691533841, 0.996181357106744};
    for (int i = 1; i < 4; ++i)
    {
      double val = update->compute(keff[i], lambda[i-1], J);
      printf("%16.12f %16.12f %16.12f \n", keff[i], val, ref[i]);
      TEST(soft_equiv(val, ref[i]));
    }
  }
  ManagerERME::finalize();
  return 0;
}

//----------------------------------------------------------------------------//
int test_AnghelUpdateINVLIN(int argc, char *argv[])
{
  {
    build_manager(argc, argv);
    build_update("anghel", AnghelUpdate::INVLIN);
    double ref[4] = {keff[0], keff[1], 0.996165591538771, 0.996181359957469};
    for (int i = 1; i < 4; ++i)
    {
      double val = update->compute(keff[i], lambda[i-1], J);
      printf("%16.12f %16.12f %16.12f \n", keff[i], val, ref[i]);
      TEST(soft_equiv(val, ref[i]));
    }
  }
  ManagerERME::finalize();
  return 0;
}

//----------------------------------------------------------------------------//
int test_AnghelUpdateINVINV(int argc, char *argv[])
{
  {
    build_manager(argc, argv);
    build_update("anghel", AnghelUpdate::INVINV);
    double ref[4] = {keff[0], keff[1], 0.996174029746272, 0.996181388727932};
    double val;
    for (int i = 1; i < 4; ++i)
    {
      val = update->compute(keff[i], lambda[i-1], J);
      printf("%16.12f %16.12f %16.12f \n", keff[i], val, ref[i]);
      TEST(soft_equiv(val, ref[i]));
    }
    printf("%16.12e \n", val - keff[4]);
  }
  ManagerERME::finalize();
  return 0;
}

//----------------------------------------------------------------------------//
//              end of file test_EigenvalueUpdate.cc
//----------------------------------------------------------------------------//
