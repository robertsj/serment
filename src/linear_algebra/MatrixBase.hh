//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   MatrixBase.hh
 *  @brief  MatrixBase
 *  @note   Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef linear_algebra_MATRIXBASE_HH_
#define linear_algebra_MATRIXBASE_HH_

#include "Vector.hh"
#include "comm/Comm.hh"
#include "utilities/DBC.hh"
#include "utilities/SP.hh"
#include "petsc.h"
#include <vector>
#include <string>

namespace linear_algebra
{

class MatrixBase
{

public:

  //--------------------------------------------------------------------------//
  // ENUMERATIONS
  //--------------------------------------------------------------------------//

  enum matrix_display_types
  {
    STDOUT, ASCII, BINARY, END_MATRIX_DISPLAY_TYPES
  };

  enum matrix_assembly_types
  {
    FINAL, FLUSH, END_MATRIX_ASSEMBLY_TYPES
  };

  enum matrix_insert_mode
  {
    INSERT, ADD, END_MATRIX_INSERT_MODES
  };

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<MatrixBase>  SP_matrix;
  typedef unsigned int                      size_type;
  typedef std::vector<int>                  vec_int;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /// Default constructor, yielding a safe, self-destructing wrapper for Mat
  MatrixBase();

  /**
   *  @brief Constructor
   *  @param    m     local number of rows
   *  @param    n     local number of columns
   */
  MatrixBase(const size_type m,
             const size_type n);

  /// Virtual destructor
  virtual ~MatrixBase() = 0;

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// Set the global sizes and local row and column bounds
  void set_sizes_and_bounds();

  /**
   *  @brief Insert values
   *
   *  @param number_rows    Number of column indices
   *  @param rows           Indices of global rows
   *  @param number_columns Number of column indices
   *  @param columns        Indices of global columns
   *  @param values         Logically 2D array of values to insert
   */
  virtual void insert_values(const size_type  number_rows,
                             const int       *rows,
                             const size_type  number_columns,
                             const int       *columns,
                             const double    *values,
                             const int        mode = INSERT)
  {
    THROW("NOT IMPLEMENTED");
  }

  /// Copy in values from another matrix that must be no larger than this one
  virtual void insert_values(SP_matrix m_in, const int mode = INSERT_VALUES)
  {
    THROW("NOT IMPLEMENTED");
  }

  /// Assemble the matrix.
  virtual void assemble(const int mode = MAT_FINAL_ASSEMBLY);

  //@{
  ///  Matrix-vector multiplication and its transpose
  virtual void multiply(Vector &v_in, Vector &v_out) = 0;
  virtual void multiply_transpose(Vector &v_in, Vector &v_out) = 0;
  //@}

  virtual SP_matrix multiply(SP_matrix m_in) = 0;

  //@{
  ///  Getters
  Mat A();
  size_type number_global_rows() const;
  size_type number_global_columns() const;
  size_type number_local_rows() const;
  size_type number_local_columns() const;
  int lower_bound() const;
  int upper_bound() const;
  bool is_assembled() const;
  //@}

  /// Set the PETSc matrix
  void set_A(Mat A);

  /**
   *  @brief Display the matrix to screen (or to output)
   *  @param  output    Flag indicating (stdout=0, ascii=1, binary=2)
   *  @param  name      File name for ascii or binary file
   */
  virtual void display(const int          output = 0,
                       const std::string &name   = "matrix.out") const;

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// PETSc matrix
  Mat d_A;
  /// Number of matrix rows
  size_type d_number_global_rows;
  /// Number of matrix columns
  size_type d_number_global_columns;
  /// Number of matrix rows on this process
  size_type d_number_local_rows;
  /// Number of matrix columns on this process
  size_type d_number_local_columns;
  /// Starting (global) index of rows on this process
  size_type d_lower_bound;
  /// Bounding index from above for rows on this process
  size_type d_upper_bound;
  /// Flag for whether the matrix is assembled
  bool d_is_assembled;

};

} // end namespace linear_algebra

#endif // linear_algebra_MATRIXBASE_HH_

//----------------------------------------------------------------------------//
//              end of file MatrixBase.hh
//----------------------------------------------------------------------------//
