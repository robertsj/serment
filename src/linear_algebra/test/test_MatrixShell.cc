//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_MatrixShell.cc
 *  @brief Test of MatrixShell class.
 *  @note  Copyright (C) 2012 Jeremy Roberts.
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST         \
        FUNC(test_MatrixShell)

#include "utilities/TestDriver.hh"
#include "linear_algebra/MatrixShell.hh"
#include "linear_algebra/LinearAlgebraSetup.hh"
#include <iostream>
#include "matrix_shell_fixture.hh"

using namespace serment_comm;
using namespace linear_algebra;
using namespace detran_test;
using std::cout;
using std::endl;

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
int test_MatrixShell(int argc, char *argv[])
{
  // Get size and rank
  int size = Comm::size();
  int rank = Comm::rank();

  // Create MatrixShell
  MatrixShell::size_type n = 5;
  MyMatrixShell A(n, n);
  cout << "n = " << " size = " << size << " ngr = " << A.number_global_rows() << endl;
  TEST(A.number_local_rows() == n);
  TEST(A.number_local_columns() == n);
  TEST(A.number_global_rows() == n * size);
  TEST(A.number_global_columns() == n * size);

  // Create two vectors
  Vector X(n);
  Vector Y(n);

  for (int i = X.lower_bound(); i < X.upper_bound(); i++)
  {
    double value = i*i;
    X.insert_values(1, &i, &value);
  }
  X.assemble();
  Y.assemble();

  // Multiply
  A.multiply(X, Y);

  // Test the vector output.
  double ref[] = {1, 2, -size * 5 + 2};
  for (int i = 0; i < n; i++)
  {
    double ref = 2.0;
    if (i == 0 and rank == 0) ref = 1.0;
    if (i == n - 1 and rank == size - 1) ref = 2 - 25 * size * size;
    TEST(detran_utilities::soft_equiv(Y[i], ref));
  }

  A.display();
  X.display();
  Y.display();

  Y[0] = 2.0;
  Y.display();
  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_MatrixShell.cc
//----------------------------------------------------------------------------//
