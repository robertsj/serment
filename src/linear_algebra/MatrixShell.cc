//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MatrixShell.cc
 * \brief  MatrixShell member definitions
 * \author Jeremy Roberts
 * \date   Aug 20, 2012
 */
//---------------------------------------------------------------------------//

#include "MatrixShell.hh"

namespace linear_algebra
{

MatrixShell::MatrixShell(const size_type m,
                         const size_type n,
                         void* context)
  : MatrixBase(m, n)
{
  // No preconditions

  // Set matrix type and context
  PetscErrorCode ierr;
  ierr = MatSetType(d_A, MATSHELL);
  ierr = MatShellSetContext(d_A, context);

  // Set the local/global size along with row bounds.
  set_sizes_and_bounds();

  // There is no assembly for a shell.
  d_is_assembled = true;

  // Set the matrix shell operations.
  ierr = set_operation();

  // Postconditions
  Ensure(!ierr);
}

} // end namespace linear_algebra

//---------------------------------------------------------------------------//
//              end of file MatrixShell.cc
//---------------------------------------------------------------------------//
