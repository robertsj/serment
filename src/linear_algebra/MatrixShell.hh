//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MatrixShell.hh
 * \brief  MatrixShell class definition.
 * \author Jeremy Roberts
 * \date   Aug 20, 2012
 */
//---------------------------------------------------------------------------//

#ifndef MATRIXSHELL_HH_
#define MATRIXSHELL_HH_

// Serment Linear Algebra
#include "MatrixBase.hh"

// Detran Utilities
#include "DBC.hh"

namespace linear_algebra
{

/*!
 *  \class MatrixShell
 *  \brief Shell matrix class for user storage and operator scheme
 *
 *  To use this class, a user must inherit it and implement the
 *  relevant shell methods, e.g. shell_multiply.  Unfortunately,
 *  inside these methods the user must work entirely with the
 *  PETSc API.
 */
class MatrixShell: public MatrixBase
{

public:

  MatrixShell(const size_type m,
              const size_type n,
              void* context);

  // No inserting in a shell.
  virtual void insert_values(const size_type number_rows,
                             const int *rows,
                             const size_type number_columns,
                             const int *columns,
                             const double *values)
  {
    THROW("Can't insert into a shell matrix.");
  }

  /// Assemble the matrix.
  virtual void assemble()
  {
    THROW("Can't assemble a shell matrix.");
  }

  //---------------------------------------------------------------------------//
  // SHELL MATRIX OPERATIONS
  //---------------------------------------------------------------------------//

  /*!
   *  \brief Matrix-vector multiplication
   *  \param x  Input vector
   *  \param y  Output vector
   */
  virtual PetscErrorCode shell_multiply(Vec x, Vec y) = 0;

  /*!
   *  \brief Matrix-vector multiplication using matrix transpose.
   *  \param x  Input vector
   *  \param y  Output vector
   */
  virtual PetscErrorCode shell_multiply_transpose(Vec x, Vec y) = 0;


protected:

  /// This must be called
  PetscErrorCode set_operation();

};

//---------------------------------------------------------------------------//
// EXTERNAL WRAPPER FUNCTION AND SETTER
//---------------------------------------------------------------------------//

PetscErrorCode multiply_wrapper(Mat A, Vec x, Vec y)
{

  // Get the context and cast as MatrixShell
  PetscErrorCode ierr;
  void *context;
  ierr = MatShellGetContext(A, &context);
  Assert(!ierr);

  MatrixShell *foo = (MatrixShell*) context;
  // Call the actual apply operator.
  return foo->shell_multiply(x, y);
}

PetscErrorCode MatrixShell::set_operation()
{
  return MatShellSetOperation(d_A, MATOP_MULT,
                              (void(*)(void)) multiply_wrapper);
}

} // end namespace linear_algebra

#endif // MATRIXSHELL_HH_ 

//---------------------------------------------------------------------------//
//              end of file MatrixShell.hh
//---------------------------------------------------------------------------//
