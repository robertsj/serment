//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MatrixBase.cc
 * \brief  MatrixBase 
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 */
//---------------------------------------------------------------------------//

// Linear Algebra
#include "MatrixBase.hh"

#include <iostream>

namespace linear_algebra
{

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
  // Preconditions
  Require(m > 0);
  Require(n > 0);

  PetscErrorCode ierr;
  ierr = MatCreate(PETSC_COMM_WORLD, &d_A);
  ierr = MatSetSizes(d_A, m, n, PETSC_DETERMINE, PETSC_DETERMINE);

  // Postconditions
  Ensure(!ierr);
}

MatrixBase::~MatrixBase()
{
  MatDestroy(&d_A);
}

void MatrixBase::assemble()
{
  if (!d_is_assembled)
  {
    MatAssemblyBegin(d_A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(d_A, MAT_FINAL_ASSEMBLY);
    d_is_assembled = true;
  }
}

// Default implementation.  Note, shell matrices have to do something else!
void MatrixBase::display(const int output, const std::string name) const
{
  if (output == STDOUT)
  {
    MatView(d_A, PETSC_VIEWER_STDOUT_WORLD);
    return;
  }
  PetscViewer viewer;
  if (output == BINARY)
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, name.c_str(), FILE_MODE_WRITE, &viewer);
  else if (output == ASCII)
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, name.c_str(), &viewer);
  else if (output == MATLAB)
    PetscViewerBinaryMatlabOpen(PETSC_COMM_WORLD, name.c_str(), &viewer);
  else
    THROW("Invalid output switch for Matrix display");
  MatView(d_A, viewer);
  PetscViewerDestroy(&viewer);
}

//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//

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
  int ngr, ngc;
  ierr = MatGetSize(d_A, &ngr, &ngc);
  Assert(ngr > 0);
  Assert(ngc > 0);
  d_number_global_rows = ngr;
  d_number_global_columns = ngc;

  // Postconditions
  Ensure(!ierr);
  Ensure(d_upper_bound - d_lower_bound == d_number_local_rows);
}

} // end namespace linear_algebra

//---------------------------------------------------------------------------//
//              end of file MatrixBase.cc
//---------------------------------------------------------------------------//
