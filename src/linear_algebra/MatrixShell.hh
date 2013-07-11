//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   MatrixShell.hh
 *  @brief  MatrixShell class definition.
 *  @note   Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef linear_algebra_MATRIXSHELL_HH_
#define linear_algebra_MATRIXSHELL_HH_

#include "MatrixBase.hh"

namespace linear_algebra
{

// Forward declare the wrapper
PetscErrorCode multiply_wrapper(Mat A, Vec x, Vec y);

/**
 *  @class MatrixShell
 *  @brief Shell matrix class for user storage and operator scheme
 *
 *  To use this class, a user must inherit it and implement the
 *  relevant shell methods, e.g. shell_multiply.  Unfortunately,
 *  inside these methods the user must work entirely with the
 *  PETSc API.
 */
class MatrixShell: public MatrixBase
{

public:

  /// Constructor with void pointer to the user context
  MatrixShell(const size_type m,
              const size_type n,
              void* context);

  // No inserting in a shell.
  void insert_values(const size_type  number_rows,
                     const int       *rows,
                     const size_type  number_columns,
                     const int       *columns,
                     const double    *values)
  {
    THROW("Can't insert into a shell matrix.");
  }

  /// Assemble the matrix.
  void assemble()
  {
    /* ... */
  }

  //@{
  ///  Matrix-vector multiplication and its transpose must be implemented
  virtual void multiply(Vector &v_in, Vector &v_out) = 0;
  virtual void multiply_transpose(Vector &v_in, Vector &v_out) = 0;
  //@}

  SP_matrix multiply(SP_matrix m_in)
  {
    THROW("NOT IMPLEMENTED");
    return SP_matrix(NULL);
  }

  /**
   *  @brief Display the matrix to screen (or to output)
   *  @param  output    Flag indicating (stdout=0, ascii=1, binary=2)
   *  @param  name      File name for ascii or binary file
   */
  virtual void display(const int          output = 0,
                       const std::string &name   = "matrix.out") const;

};

//----------------------------------------------------------------------------//
// WRAPPERS FOR MULTIPLICATIONS
//----------------------------------------------------------------------------//

inline PetscErrorCode shell_multiply_wrapper(Mat A, Vec x, Vec y)
{
  // get the context and cast
  PetscErrorCode ierr;
  void *context;
  ierr = MatShellGetContext(A, &context);
  Assert(!ierr);
  MatrixShell* foo = (MatrixShell*) context;
  // wrap the petsc vectors
  Vector X(x), Y(y);
  // call the actual apply operator.
  foo->multiply(X, Y);
  return ierr;
}

inline PetscErrorCode shell_multiply_transpose_wrapper(Mat A, Vec x, Vec y)
{
  // get the context and cast
  PetscErrorCode ierr;
  void *context;
  ierr = MatShellGetContext(A, &context);
  Assert(!ierr);
  MatrixShell* foo = (MatrixShell*) context;
  // wrap the petsc vectors
  Vector X(x), Y(y);
  // call the actual apply operator.
  foo->multiply_transpose(X, Y);
  return ierr;
}



} // end namespace linear_algebra

#endif // linear_algebra_MATRIXSHELL_HH_

//----------------------------------------------------------------------------//
//              end of file MatrixShell.hh
//----------------------------------------------------------------------------//
