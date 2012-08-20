//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Matrix.cc
 * \author robertsj
 * \date   Aug 20, 2012
 * \brief  Matrix member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#include "Matrix.hh"

namespace linear_algebra
{

Matrix::Matrix(const size_type m,
               const size_type n,
               const vec_int &number_nonzeros,
               const vec_int &number_nonzeros_off)
  : MatrixBase(m, n)
{
  // Preconditions
  Require(number_nonzeros.size() > 0);
  Require(number_nonzeros_off.size() > 0);

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

  // Postconditions
  Ensure(!ierr);
}


} // end namespace linear_algebra

//---------------------------------------------------------------------------//
//              end of file Matrix.cc
//---------------------------------------------------------------------------//


