//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   Matrix.cc
 *  @brief  Matrix member definitions
 *  @note   Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "Matrix.hh"
#include <iostream>

namespace linear_algebra
{

//----------------------------------------------------------------------------//
Matrix::Matrix(const size_type  m,
               const size_type  n,
               const vec_int   &number_nonzeros,
               const vec_int   &number_nonzeros_off)
  : MatrixBase(m, n)
{
  Require(number_nonzeros.size() > 0);
  Require(number_nonzeros_off.size() > 0);

  preallocate(number_nonzeros, number_nonzeros_off);
}

//----------------------------------------------------------------------------//
Matrix::Matrix(const size_type m,
               const size_type n)
  : MatrixBase(m, n)
{
  /* ... */
}

//----------------------------------------------------------------------------//
void Matrix::insert_values(const size_type   number_rows,
                           const int        *rows,
                           const size_type   number_columns,
                           const int        *columns,
                           const double     *values,
                           const int         mode)
{
  Require(number_rows > 0);
  Require(number_columns > 0);

  InsertMode petsc_mode = INSERT_VALUES;
  if (mode == ADD) petsc_mode = ADD_VALUES;
  d_is_assembled = false; // changing values = no longer ready
  PetscErrorCode ierr;
  ierr = MatSetValues(d_A, number_rows, rows, number_columns, columns,
                      values, petsc_mode);

  Ensure(!ierr);
}

//----------------------------------------------------------------------------//
void Matrix::insert_values(SP_matrix M_in, const int mode)
{
  Insist(dynamic_cast<Matrix*>(M_in.bp()), "Values inserted into Matrix only");
  Require(M_in->number_local_columns() <= number_local_columns());
  Require(M_in->number_local_rows()    <= number_local_rows());

  Matrix &M = *dynamic_cast<Matrix*>(M_in.bp());

  int ncols = 0;
  const int* cols = NULL;
  const double *vals = NULL;
  for (int i = M.lower_bound(); i < M.upper_bound(); ++i)
  {
    MatGetRow(M.A(), i, &ncols, &cols, &vals);
    int rows[] = {i};
    insert_values(1, rows, ncols, cols, vals, mode);
    MatRestoreRow(M.A(), i, &ncols, &cols, &vals);
  }
}

//----------------------------------------------------------------------------//
void Matrix::preallocate(const vec_int &nnz,
                         const vec_int &nnz_od)
{
  Require(nnz.size() > 0);
  Require(nnz.size() == nnz_od.size());

  // Create appropriate matrix type.
  PetscErrorCode ierr;
  int size; // \todo Why is MPI called directly?
  //ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);

  size = serment_comm::Comm::size();

  if (size > 1)
  {
    ierr = MatSetType(d_A, MATMPIAIJ);
    ierr = MatMPIAIJSetPreallocation(d_A,
                                     PETSC_NULL, &nnz[0],
                                     PETSC_NULL, &nnz_od[0]);
  }
  else
  {
    // All nonzeros get lumped into one vector.
    vec_int total_nnz(nnz);
    for (int i = 0; i < nnz.size(); ++i)
    {
      total_nnz[i] += nnz_od[i];
      //std::cout << " i ---> " << total_nnz[i] << std::endl;
    }
    ierr = MatSetType(d_A, MATSEQAIJ);
    ierr = MatSeqAIJSetPreallocation(d_A, PETSC_NULL, &total_nnz[0]);
  }

  ierr = MatSetOption(d_A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

  // Set the local/global size along with row bounds.
  set_sizes_and_bounds();

  Ensure(!ierr);
}

//----------------------------------------------------------------------------//
Matrix::SP_matrix Matrix::multiply(SP_matrix m_in)
{
  SP_matrix m_out(new Matrix);
  Mat m_out_A;
  MatMatMult(d_A, m_in->A(), MAT_INITIAL_MATRIX, 1.0, &m_out_A);
  m_out->set_A(m_out_A);
  return m_out;
}

} // end namespace linear_algebra

//----------------------------------------------------------------------------//
//              end of file Matrix.cc
//----------------------------------------------------------------------------//


