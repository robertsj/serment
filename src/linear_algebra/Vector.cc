//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Vector.cc
 * \brief  Vector 
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 */
//---------------------------------------------------------------------------//

#include "Vector.hh"

namespace linear_algebra
{

Vector::Vector(const size_type m)
  : d_local_size(m)
  , d_is_assembled(false)
{
  // Preconditions
  Require(m > 0);

  // Create the vector
  PetscErrorCode ierr;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, m, PETSC_DETERMINE, &d_V);

  // Get the global size.
  int gs;
  ierr = VecGetSize(d_V, &gs);
  Assert(gs > 0);
  d_global_size = gs;

  // Get the global ranges.
  int lb, ub;
  ierr = VecGetOwnershipRange(d_V, &lb, &ub);
  Assert(lb >= 0 and ub > 0);
  d_lower_bound = lb;
  d_upper_bound = ub;

  // Get acces to and then "restore" the array, but not actually.
  // PETSc requires the restore call to ensure
  // safe coding, but by passing null, we get to keep it.
  // We'll code safely...
  ierr = VecGetArray(d_V, &d_array);
  ierr = VecRestoreArray(d_V, PETSC_NULL);

  // Postconditions
  Ensure(!ierr);
  Ensure(d_upper_bound - d_lower_bound == d_local_size);
}

Vector::~Vector()
{
  VecDestroy(&d_V);
}

void Vector::assemble()
{
  if (!d_is_assembled)
  {
    VecAssemblyBegin(d_V);
    VecAssemblyEnd(d_V);
    d_is_assembled = true;
  }
}

void Vector::display() const
{
  VecView(d_V, PETSC_VIEWER_STDOUT_SELF);
}

} // end namespace linear_algebra

//---------------------------------------------------------------------------//
//              end of file Vector.cc
//---------------------------------------------------------------------------//
