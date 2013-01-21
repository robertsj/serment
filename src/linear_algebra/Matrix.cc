//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Matrix.cc
 *  @author robertsj
 *  @date   Aug 20, 2012
 *  @brief  Matrix member definitions.
 */
//---------------------------------------------------------------------------//

#include "Matrix.hh"
#include <iostream>
namespace linear_algebra
{

//---------------------------------------------------------------------------//
Matrix::Matrix(const size_type m,
               const size_type n,
               const vec_int &number_nonzeros,
               const vec_int &number_nonzeros_off)
  : MatrixBase(m, n)
{
  // Preconditions
  Require(number_nonzeros.size() > 0);
  Require(number_nonzeros_off.size() > 0);

  // Create appropriate matrix type.
  PetscErrorCode ierr;
  int size;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);
  if (size > 1)
  {
    ierr = MatSetType(d_A, MATMPIAIJ);
    ierr = MatMPIAIJSetPreallocation(d_A,
                                     PETSC_NULL, &number_nonzeros[0],
                                     PETSC_NULL, &number_nonzeros_off[0]);
  }
  else
  {
    ierr = MatSetType(d_A, MATSEQAIJ);
    ierr = MatSeqAIJSetPreallocation(d_A, PETSC_NULL, &number_nonzeros[0]);
  }

  // Set the local/global size along with row bounds.
  set_sizes_and_bounds();

  // Postconditions
  Ensure(!ierr);
}

//---------------------------------------------------------------------------//
Matrix::Matrix(const size_type m,
               const size_type n)
  : MatrixBase(m, n)
{
  /* ... */
}

//---------------------------------------------------------------------------//
void Matrix::
preallocate(const vec_int &number_nonzeros,
            const vec_int &number_nonzeros_off)
{
  // Preconditions
  Require(number_nonzeros.size() > 0);
  Require(number_nonzeros.size() == number_nonzeros_off.size());

  // Create appropriate matrix type.
  PetscErrorCode ierr;
  int size;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);
  if (size > 1)
  {
    ierr = MatSetType(d_A, MATMPIAIJ);
    ierr = MatMPIAIJSetPreallocation(d_A,
                                     PETSC_NULL, &number_nonzeros[0],
                                     PETSC_NULL, &number_nonzeros_off[0]);
  }
  else
  {
    // All nonzeros get lumped into one vector.
    vec_int total_number_nonzeros(number_nonzeros);
    for (size_t i = 0; i < number_nonzeros.size(); ++i)
      total_number_nonzeros[i] += number_nonzeros_off[i];

    ierr = MatSetType(d_A, MATSEQAIJ);
    ierr = MatSeqAIJSetPreallocation(d_A, PETSC_NULL, &total_number_nonzeros[0]);
  }

  // Set the local/global size along with row bounds.
  set_sizes_and_bounds();

  // Postconditions
  Ensure(!ierr);
}



} // end namespace linear_algebra

//---------------------------------------------------------------------------//
//              end of file Matrix.cc
//---------------------------------------------------------------------------//


