//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   Matrix.hh
 *  @brief  Matrix class definition
 *  @note   Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef linear_algebra_MATRIX_HH_
#define linear_algebra_MATRIX_HH_

#include "MatrixBase.hh"

namespace linear_algebra
{

/// Lightweight wrapper around a PETSc compressed sparse row matrix
class Matrix: public MatrixBase
{

public:

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param    m       local number of rows
   *  @param    n       local number of columns
   *  @param    nnz     number of nonzeros per row in the local columns
   *  @param    nnz_od  number of nonzeros per row out of the local columns
   */
  Matrix(const size_type  m,
         const size_type  n,
         const vec_int   &nnz,
         const vec_int   &nnz_od = vec_int(0));

  /// Constructor with deferred preallocation
  Matrix(const size_type m,
         const size_type n);


  /// Default constructor
  Matrix(){}

  /// Insert values explicitly
  void insert_values(const size_type  number_rows,
                     const int       *rows,
                     const size_type  number_columns,
                     const int       *columns,
                     const double    *values,
                     const int        mode = INSERT);

  /// Insert values by copying a possibly smaller matrix
  void insert_values(SP_matrix m_in, const int mode = INSERT);

  /// Preallocate the matrix with the number of nonzeros
  void preallocate(const vec_int &nnz,
                   const vec_int &nnz_od = vec_int(0));

  //@{
  ///  Use the base implementation for matrix-vector operations
  void multiply(Vector &v_in, Vector &v_out)
  {
    MatrixBase::multiply(v_in, v_out);
  }
  void multiply_transpose(Vector &v_in, Vector &v_out)
  {
    MatrixBase::multiply_transpose(v_in, v_out);
  }
  //@}

  ///  Matrix-matrix multiplication
  SP_matrix multiply(SP_matrix m_in);

};

} // end namespace linear_algebra

#endif // linear_algebra_MATRIX_HH_

//----------------------------------------------------------------------------//
//              end of file Matrix.hh
//----------------------------------------------------------------------------//
