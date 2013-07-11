//----------------------------------*-C++-*-----------------------------------//
/**
 * @file  test_Matrix.cc
 * @brief Test of Matrix class.
 * @note  Copyright (C) 2012 Jeremy Roberts.
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                   \
        FUNC(test_Matrix)           \
        FUNC(test_Matrix_matmatmult)

#include "utilities/TestDriver.hh"
#include "linear_algebra/Matrix.hh"
#include "linear_algebra/LinearAlgebraSetup.hh"
#include <iostream>

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
int test_Matrix(int argc, char *argv[])
{
  // Get size and rank
  int size = Comm::size();
  int rank = Comm::rank();

  // Create Matrix
  Matrix::size_type n = 5;
  Matrix::vec_int dnnz(5, 3);  // At most 3 elements in local block
  Matrix::vec_int odnnz(5, 1); // At most 1 element outside local block
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
int test_Matrix_matmatmult(int argc, char *argv[])
{
  int m_A = 5 * Comm::size();
  int n_A = 5 * Comm::size();
  int m_B = m_A;
  int n_B = n_A ;
  Matrix::vec_int nnz_d_A(m_A, 3);
  Matrix::vec_int nnz_o_A(m_A, 0);
  Matrix::vec_int nnz_d_B(m_B, 3);
  Matrix::vec_int nnz_o_B(m_B, 0);
  Matrix::SP_matrix A(new Matrix(m_A, n_A, nnz_d_A, nnz_o_A));
  Matrix::SP_matrix B(new Matrix(m_B, n_B, nnz_d_B, nnz_o_B));

  for (int i = A->lower_bound(); i < A->upper_bound(); ++i)
  {
    if (i == 0)
    {
      double v[] = {2.0, -1.0};
      int r[] = {0};
      int c[] = {0, 1};
      A->insert_values(1, r, 2, c, v);
      B->insert_values(2, c, 1, r, v);
    }
    else if (i == m_A - 1)
    {
      double v[] = {-1.0, 2.0};
      int r[] = {i};
      int c[] = {i - 1, i};
      A->insert_values(1, r, 2, c, v);
      B->insert_values(2, c, 1, r, v);
    }
    else
    {
      double v[] = {-1.0, 2.0, -1.0};
      int r[] = {i};
      int c[] = {i - 1, i, i + 1};
      A->insert_values(1, r, 3, c, v);
      B->insert_values(3, c, 1, r, v);
    }
  }
  A->assemble();
  A->display();
  B->assemble();
  B->display();

  Matrix::SP_matrix C = B->multiply(A);
  return 0;
}
//----------------------------------------------------------------------------//
//              end of test_Matrix.cc
//----------------------------------------------------------------------------//
