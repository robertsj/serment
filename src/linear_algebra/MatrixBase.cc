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
  , d_is_assembled(false)
{
  using std::cout;
  using std::endl;
  // Preconditions
  Require(m > 0);
  Require(n > 0);

  PetscErrorCode ierr;
  ierr = MatCreate(PETSC_COMM_WORLD, &d_A);
  ierr = MatSetSizes(d_A, m, n, PETSC_DETERMINE, PETSC_DETERMINE);

  // Postconditions
  Require(!ierr);
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
void MatrixBase::display() const
{
  MatView(d_A, PETSC_VIEWER_STDOUT_WORLD);
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
  Require(!ierr);
  Ensure(d_upper_bound - d_lower_bound == d_number_local_rows);
}

} // end namespace linear_algebra

//---------------------------------------------------------------------------//
//              end of file MatrixBase.cc
//---------------------------------------------------------------------------//
