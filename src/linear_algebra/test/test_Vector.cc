//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   test_Vector.cc
 *  @brief  Test of Vector class.
 *  @note   Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                         \
        FUNC(test_Vector)                 \
        FUNC(test_Vector_temporary)       \
        FUNC(test_Vector_collect_to_root)

#include "utilities/TestDriver.hh"
#include "Vector.hh"
#include "LinearAlgebraSetup.hh"
#include <iostream>

using namespace serment_comm;
using namespace linear_algebra;
using namespace detran_test;
using namespace detran_utilities;
using std::cout;
using std::endl;
#define COUT(c) cout << c << endl;

int main(int argc, char *argv[])
{
  linear_algebra::initialize(argc, argv, true);
  RUN(argc, argv);
  linear_algebra::finalize();
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

// Test of basic public interface
int test_Vector(int argc, char *argv[])
{
  // Create vector
  Vector X(10);
  TEST(X.local_size() == 10);

  double value = 1.0;
  int count = 1;
  for (int i = X.lower_bound(); i < X.upper_bound(); i++)
    X.insert_values(count, &i, &value);
  X.assemble();

  value = X.dot(X);
  TEST(detran_utilities::soft_equiv(value, 1.0*X.global_size()));

  X.scale(2.0);
  for (int i = 0; i < X.local_size(); i++)
    TEST(detran_utilities::soft_equiv(X[i], 2.0));

  Vector V(X);
  for (int i = 0; i < V.local_size(); i++)
    TEST(detran_utilities::soft_equiv(V[i], 2.0));

  return 0;
}

// Test of temporary constructors
int test_Vector_temporary(int argc, char *argv[])
{
  Vector X(5, 0.0);
  for (int i = 0; i < X.local_size(); ++i)
  {
    X[i] = Comm::rank() * 5.0 + i;
  }
  double sum_ref = (X.global_size()-1) * X.global_size() /2;
  double sum_X   = X.norm(X.L1);
  TEST(soft_equiv(sum_X, sum_ref));

  // test temporary from PETSc Vec
  {
    Vector temp_from_Vec(X.V());
    double sum_temp = temp_from_Vec.norm(X.L1);
    TEST(soft_equiv(sum_temp, sum_ref));
  }

  // test temporary by inserting larger vec into smaller vec
  {
    int m = 5;
    if (Comm::rank() == Comm::size() - 1) m = 3;
    Vector small(X, m);
    TEST(small.global_size() == X.global_size() - 2);
    double sum_small = small.norm(X.L1);
    sum_ref -= (2*X.global_size()-3);
    TEST(soft_equiv(sum_small, sum_ref));
  }

  return 0;
}

// Test of collect to a root
int test_Vector_collect_to_root(int argc, char *argv[])
{
  using serment_comm::Comm;
  Vector V(1, 0.0);

  // test ranges
  Vector::vec_int ranges = V.ranges();
  for (int i = 0; i < ranges.size(); ++i)
  {
    TEST(ranges[i] == i);
  }

  // fill with local values
  for (int i = 0; i < V.local_size(); ++i)
  {
    V[i] = Comm::rank();
  }

  // collect on last process
  int root = Comm::last();
  Vector::SP_vector V_seq = V.collect_on_root(root);
  if (Comm::rank() == root)
  {
    TEST(V_seq->is_sequential());
    TEST(V_seq->local_size() == V_seq->global_size());
    for (int i = 0; i < V_seq->global_size(); ++i)
    {
      TEST(soft_equiv((*V_seq)[i], (double)i));
    }
  }

  Comm::global_barrier();
  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_Vector.cc
//----------------------------------------------------------------------------//
