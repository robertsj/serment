//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   MatrixBase.cc
 *  @brief  MatrixBase
 *  @note   Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "MatrixBase.hh"
#include <iostream>

namespace linear_algebra
{

//----------------------------------------------------------------------------//
MatrixBase::MatrixBase()
  : d_A(PETSC_NULL)
  , d_number_local_rows(0)
  , d_number_local_columns(0)
  , d_number_global_rows(0)
  , d_number_global_columns(0)
  , d_lower_bound(0)
  , d_upper_bound(0)
  , d_is_assembled(false)
{
  /* ... */
}

//----------------------------------------------------------------------------//
MatrixBase::MatrixBase(const size_type m,
                       const size_type n)
  : d_number_local_rows(m)
  , d_number_local_columns(n)
  , d_number_global_rows(0)
  , d_number_global_columns(0)
  , d_lower_bound(0)
  , d_upper_bound(0)
  , d_is_assembled(false)
{
  Require(m > 0);
  Require(n > 0);

  PetscErrorCode ierr;
  ierr = MatCreate(PETSC_COMM_WORLD, &d_A);
  ierr = MatSetSizes(d_A, m, n, PETSC_DETERMINE, PETSC_DETERMINE);

  Ensure(!ierr);
}

//----------------------------------------------------------------------------//
MatrixBase::~MatrixBase()
{
  if (d_A != PETSC_NULL) MatDestroy(&d_A);
}

//----------------------------------------------------------------------------//
void MatrixBase::assemble(const int mode)
{
  if (!d_is_assembled)
  {
    if (mode == FINAL)
    {
      MatAssemblyBegin(d_A, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(d_A, MAT_FINAL_ASSEMBLY);
    }
    else
    {
      MatAssemblyBegin(d_A, MAT_FLUSH_ASSEMBLY);
      MatAssemblyEnd(d_A, MAT_FLUSH_ASSEMBLY);
    }
    d_is_assembled = true;
  }
}

//----------------------------------------------------------------------------//
void MatrixBase::multiply(Vector &v_in, Vector &v_out)
{
  Require(d_is_assembled);
  Require(d_number_global_columns = v_in.global_size());
  Require(d_number_global_rows = v_out.global_size());

  PetscErrorCode ierr = MatMult(d_A, v_in.V(), v_out.V());

  Ensure(!ierr);
}

//----------------------------------------------------------------------------//
void MatrixBase::multiply_transpose(Vector &v_in, Vector &v_out)
{
  Require(d_is_assembled);
  Require(d_number_global_columns = v_out.global_size());
  Require(d_number_global_rows = v_in.global_size());

  PetscErrorCode ierr = MatMultTranspose(d_A, v_in.V(), v_out.V());

  Ensure(!ierr);
}

//----------------------------------------------------------------------------//
Mat MatrixBase::A()
{
  return d_A;
}

//----------------------------------------------------------------------------//
MatrixBase::size_type MatrixBase::number_global_rows() const
{
  return d_number_global_rows;
}

//----------------------------------------------------------------------------//
MatrixBase::size_type MatrixBase::number_global_columns() const
{
  return d_number_global_columns;
}

//----------------------------------------------------------------------------//
MatrixBase::size_type MatrixBase::number_local_rows() const
{
  return d_number_local_rows;
}

//----------------------------------------------------------------------------//
MatrixBase::size_type MatrixBase::number_local_columns() const
{
  return d_number_local_columns;
}

//----------------------------------------------------------------------------//
int MatrixBase::lower_bound() const
{
  return d_lower_bound;
}

//----------------------------------------------------------------------------//
int MatrixBase::upper_bound() const
{
  return d_upper_bound;
}

//----------------------------------------------------------------------------//
bool MatrixBase::is_assembled() const
{
  return d_is_assembled;
}

//----------------------------------------------------------------------------//
void MatrixBase::display(const int output, const std::string &name) const
{
  Require(output < Vector::END_VECTOR_DISPLAY_TYPES);
  PetscErrorCode ierr = 0;
  if (output == Vector::STDOUT)
  {
    ierr = MatView(d_A, PETSC_VIEWER_STDOUT_WORLD);
  }
  else
  {
    PetscViewer viewer;
    if (output == Vector::BINARY)
    {
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, name.c_str(),
                                   FILE_MODE_WRITE, &viewer);
    }
    else if (output == Vector::ASCII)
    {
      ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, name.c_str(), &viewer);
    }
    ierr = MatView(d_A, viewer);
    ierr = PetscViewerDestroy(&viewer);
  }
  Ensure(!ierr);
  return;
}

void MatrixBase::set_A(Mat A)
{
  d_A = A;
  set_sizes_and_bounds();
}

//----------------------------------------------------------------------------//
void MatrixBase::set_sizes_and_bounds()
{
  // Get ranges
  PetscErrorCode ierr;
  int lb, ub;

  ierr = MatGetOwnershipRange(d_A, &lb, &ub);
  Assert(lb >= 0 and ub > 0);
  d_lower_bound = lb;
  d_upper_bound = ub;

  // Get global sizes.
  int m, n;
  ierr = MatGetSize(d_A, &m, &n);
  Assert(m > 0);
  Assert(n > 0);
  d_number_global_rows = m;
  d_number_global_columns = n;
  ierr = MatGetLocalSize(d_A, &m, &n);
  Assert(m > 0);
  Assert(n > 0);
  d_number_local_rows = m;
  d_number_local_columns = n;

  Ensure(!ierr);
  Ensure(d_upper_bound - d_lower_bound == d_number_local_rows);
}

} // end namespace linear_algebra

//----------------------------------------------------------------------------//
//              end of file MatrixBase.cc
//----------------------------------------------------------------------------//
