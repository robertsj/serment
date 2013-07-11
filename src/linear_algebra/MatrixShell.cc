//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  MatrixShell.cc
 *  @brief MatrixShell member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#include "MatrixShell.hh"
#include "Matrix.hh"

namespace linear_algebra
{

//---------------------------------------------------------------------------//
MatrixShell::MatrixShell(const size_type  m,
                         const size_type  n,
                         void            *context)
  : MatrixBase(m, n)
{
  PetscErrorCode ierr;
  ierr = MatSetType(d_A, MATSHELL);
  ierr = MatShellSetContext(d_A, context);
  ierr = MatSetUp(d_A);
  set_sizes_and_bounds();

  // There is no assembly for a shell.
  d_is_assembled = true;

  ierr = MatShellSetOperation(d_A, MATOP_MULT,
                              (void(*)(void))shell_multiply_wrapper);
  ierr = MatShellSetOperation(d_A, MATOP_MULT_TRANSPOSE,
                              (void(*)(void))shell_multiply_wrapper);

  Ensure(!ierr);
}

//----------------------------------------------------------------------------//
void MatrixShell::display(const int output, const std::string &name) const
{
  Require(output < Vector::END_VECTOR_DISPLAY_TYPES);
  PetscErrorCode ierr = 0;
  Matrix M;
  Mat M_A;
  ierr = MatComputeExplicitOperator(d_A, &M_A);
  M.set_A(M_A);
  M.display(output, name);
  return;
}

} // end namespace linear_algebra

//---------------------------------------------------------------------------//
//              end of file MatrixShell.cc
//---------------------------------------------------------------------------//
