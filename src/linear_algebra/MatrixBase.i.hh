//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MatrixBase.i.hh
 * \brief  MatrixBase.i 
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 */
//---------------------------------------------------------------------------//

#ifndef MATRIXBASE_I_HH_
#define MATRIXBASE_I_HH_

namespace linear_algebra
{

// Insertion

inline void MatrixBase::insert_values(const size_type number_rows,
                                      const int *rows,
                                      const size_type number_columns,
                                      const int *columns,
                                      const double *values)
{
  // Preconditions
  Require(number_rows > 0);
  Require(number_columns > 0);

  d_is_assembled = false; // changing values = no longer ready

  PetscErrorCode ierr;
  ierr = MatSetValues(d_A, number_rows, rows, number_columns, columns,
                      values, INSERT_VALUES);

  // Postconditions
  Ensure(!ierr);
}

// Operations

inline void MatrixBase::multiply(Vector &x, Vector &y)
{
  // Preconditions
  Require(d_is_assembled);
  Require(d_number_global_columns = x.global_size());
  Require(d_number_global_rows = y.global_size());

  PetscErrorCode ierr;
  ierr = MatMult(d_A, x.V(), y.V());

  // Postconditions
  Ensure(!ierr);
}

inline void MatrixBase::multiply_transpose(Vector &x, Vector &y)
{
  // Preconditions
  Require(d_is_assembled);
  Require(d_number_global_columns = y.global_size());
  Require(d_number_global_rows = x.global_size());

  PetscErrorCode ierr;
  ierr = MatMultTranspose(d_A, x.V(), y.V());

  // Postconditions
  Ensure(!ierr);
}

} // end namespace linear_algebra

#endif // MATRIXBASE_I_HH_ 

//---------------------------------------------------------------------------//
//              end of file MatrixBase.i.hh
//---------------------------------------------------------------------------//
