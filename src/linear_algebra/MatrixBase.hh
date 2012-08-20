//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MatrixBase.hh
 * \brief  MatrixBase 
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 */
//---------------------------------------------------------------------------//

#ifndef MATRIXBASE_HH_
#define MATRIXBASE_HH_

// Linear Algebra
#include "Vector.hh"

// Detran Utilities
#include "DBC.hh"

// System
#include "petsc.h"
#include <vector>

namespace linear_algebra
{

class MatrixBase
{

public:

  typedef unsigned int size_type;
  typedef std::vector<int> vec_int;

  /*!
   *  \brief Constructor
   *  \param m
   *  \param n
   *  \param number_nonzeros
   *  \param number_nonzeros_offdiagonal
   */
  MatrixBase(const size_type m,
             const size_type n,
             const vec_int &number_nonzeros,
             const vec_int &number_nonzeros_offdiagonal = vec_int(0));

  /// Virtual destructor
  virtual ~MatrixBase();

  /*!
   *  \brief Insert values
   *  \param values         Logically 2D array of values to insert
   *  \param number_rows    Number of column indices
   *  \param rows           Indices of rows
   *  \param number_columns Number of column indices
   *  \param columns        Indices of columns
   */
  virtual void insert_values(const double *values,
                             const size_type number_rows,
                             const int *rows,
                             const size_type number_columns,
                             const int *columns);


  /// Assemble the matrix.
  virtual void assemble();

  /// \name Matrix Operations
  /// \{

  /*!
   *  \brief Matrix-vector multiplication
   *  \param x  Input vector
   *  \param y  Output vector
   */
  virtual void multiply(Vector &x, Vector &y);

  /*!
   *  \brief Matrix-vector multiplication using matrix transpose.
   *  \param x  Input vector
   *  \param y  Output vector
   */
  virtual void multiply_transpose(Vector &x, Vector &y);

  /// \}

  size_type number_global_rows() const
  {
    return d_number_global_rows;
  }

  size_type number_global_columns() const
  {
    return d_number_global_columns;
  }

  size_type number_local_rows() const
  {
    return d_number_local_rows;
  }

  size_type number_local_columns() const
  {
    return d_number_local_columns;
  }

  bool is_assembled() const
  {
    return d_is_assembled;
  }

  virtual void display() const;

private:

  /// \name Private Data
  /// \{

  /// PETSc matrix
  Mat d_A;

  const size_type d_number_global_rows;

  const size_type d_number_global_columns;

  const size_type d_number_local_rows;

  const size_type d_number_local_columns;

  bool d_is_assembled;

  /// \}


};

} // end namespace linear_algebra

// Inline members
#include "MatrixBase.i.hh"

#endif // MATRIXBASE_HH_ 

//---------------------------------------------------------------------------//
//              end of file MatrixBase.hh
//---------------------------------------------------------------------------//
