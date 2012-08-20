//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_Matrix.cc
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 * \brief  Test of Matrix class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST         \
        FUNC(test_Matrix)

// Detran test
#include "TestDriver.hh"

#include "Matrix.hh"

#include <iostream>

// Setup

using namespace linear_algebra;
using namespace detran_test;
using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_Matrix_actual();

// Test of basic public interface
int test_Matrix(int argc, char *argv[])
{
  // Initialize PETSc
  PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

  // Call actual test.
  int result = test_Matrix_actual();

  // Finalize PETSc
  PetscFinalize();

  return result;
}

// Test of basic public interface
int test_Matrix_actual()
{
  // Get size and rank
  int size, rank;
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // Create Matrix
  Matrix::size_type n = 5;
  Matrix::vec_int dnnz(5, 3);
  Matrix::vec_int odnnz(5, 0); // Nothing outside of my row-column block
  Matrix A(n, n, dnnz, odnnz);
  TEST(A.number_local_rows() == n);
  TEST(A.number_local_columns() == n);

  cout << A.number_global_rows() << " "  << n * size << endl;

  TEST(A.number_global_rows() == n * size);
  TEST(A.number_global_columns() == n * size);

  // Doing a small 1D second order difference
  for (int row = A.lower_bound(); row < A.upper_bound(); row++)
  {
    if (row == 0)
    {
      int c[] = {0, 1};
      double v[] = {-2, 1};
      A.insert_values(1, &row, 2, c, v);
    }
    else if (row == A.number_global_rows() - 1)
    {
      int c[] = {A.number_global_rows() - 2, A.number_global_rows() - 1};
      double v[] = {1, -2};
      A.insert_values(1, &row, 2, c, v);
    }
    else
    {
      int c[] = {row - 1, row, row + 1};
      double v[] = {1, -2, 1};
      A.insert_values(1, &row, 3, c, v);
    }
  }
  A.assemble();

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

//  // Test the vector output.
//  double ref[] = {1, 2, 2, 2, -23};
//  for (int i = 0; i < n; i++)
//    TEST(detran::soft_equiv(Y[i], ref[i]));

  A.display();
  X.display();
  Y.display();

  Y[0] = 2.0;
  Y.display();
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_Matrix.cc
//---------------------------------------------------------------------------//
