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

// Comm
#include "comm/Comm.hh"

// Detran Utilities
#include "DBC.hh"
#include "SP.hh"

// System
#include "petsc.h"
#include <vector>
#include <string>

namespace linear_algebra
{

class MatrixBase
{

public:

  /// Output flag options for matrix viewer
  enum MATOUT
  {
    STDOUT, ASCII, BINARY, MATLAB
  };

  typedef detran::SP<MatrixBase>  SP_matrix;
  typedef unsigned int            size_type;
  typedef std::vector<int>        vec_int;

  /*!
   *  \brief Constructor
   *  \param m local number of rows
   *  \param n local number of columns
   *  \param number_nonzeros
   *  \param number_nonzeros_offdiagonal
   */
  MatrixBase(const size_type m,
             const size_type n);

  /// Virtual destructor
  virtual ~MatrixBase() = 0;

  /*!
   *  \brief Insert values
   *
   *  \param number_rows    Number of column indices
   *  \param rows           Indices of global rows
   *  \param number_columns Number of column indices
   *  \param columns        Indices of global columns
   *  \param values         Logically 2D array of values to insert
   */
  virtual void insert_values(const size_type number_rows,
                             const int *rows,
                             const size_type number_columns,
                             const int *columns,
                             const double *values);

  /// Assemble the matrix.
  virtual void assemble();

  //---------------------------------------------------------------------------//
  // MATRIX OPERATIONS
  //---------------------------------------------------------------------------//

  /*!
   *  \brief Matrix-vector multiplication
   *  \param x  Input vector
   *  \param y  Output vector
   */
  void multiply(Vector &x, Vector &y);

  /*!
   *  \brief Matrix-vector multiplication using matrix transpose.
   *  \param x  Input vector
   *  \param y  Output vector
   */
  void multiply_transpose(Vector &x, Vector &y);


  //---------------------------------------------------------------------------//
  // GETTERS
  //---------------------------------------------------------------------------//

  /// Get PETSc Mat object.
  Mat A()
  {
    return d_A;
  }

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

  /// Return the lower bound.
  int lower_bound() const
  {
    return d_lower_bound;
  }

  /// Return the upper bound.
  int upper_bound() const
  {
    return d_upper_bound;
  }

  bool is_assembled() const
  {
    return d_is_assembled;
  }

  /*!
   *  \brief Display the matrix to screen (or to output)
   *  \param  output    Flag indicating (stdout=0, ascii=1, binary=2)
   *  \param  name      File name for ascii or binary file
   */
  virtual void display(const int output = 0,
                       const std::string name = "matrix.out") const;

protected:

  /// \name Private Data
  /// \{

  /// PETSc matrix
  Mat d_A;

  /// Number of matrix rows
  size_type d_number_global_rows;

  /// Number of matrix columns
  size_type d_number_global_columns;

  /// Number of matrix rows on this process
  const size_type d_number_local_rows;

  /// Number of matrix columns on this process
  const size_type d_number_local_columns;

  /// Starting (global) index of rows on this process
  size_type d_lower_bound;

  /// Bounding index from above for rows on this process
  size_type d_upper_bound;

  /// Flag for whether the matrix is assembled
  bool d_is_assembled;

  /// \}

  /// Set the global sizes and local row and column bounds
  void set_sizes_and_bounds();

};

} // end namespace linear_algebra

// Inline members
#include "MatrixBase.i.hh"

#endif // MATRIXBASE_HH_ 

//---------------------------------------------------------------------------//
//              end of file MatrixBase.hh
//---------------------------------------------------------------------------//
