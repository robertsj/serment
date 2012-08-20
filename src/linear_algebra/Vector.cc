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

Vector::Vector(const unsigned int m)
  : d_is_assembled(false)
  , d_global_size(m)
  , d_local_size(m)
{
  // Preconditions
  Require(m > 0);

  PetscErrorCode ierr;
  ierr = VecCreateSeq(PETSC_COMM_SELF, m, &d_V);
  Insist(!ierr, "Error creating Vec");
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
