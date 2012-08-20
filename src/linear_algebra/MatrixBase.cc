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

namespace linear_algebra
{

MatrixBase::MatrixBase(const size_type m,
                       const size_type n,
                       const vec_int &number_nonzeros,
                       const vec_int &number_nonzeros_offdiagonal)
  : d_number_global_rows(m)
  , d_number_global_columns(n)
  , d_number_local_rows(m)
  , d_number_local_columns(n)
  , d_is_assembled(false)
{
  // Preconditions
  Require(m > 0);
  Require(n > 0);
}

MatrixBase::~MatrixBase()
{
  MatDestroy(&d_A);
}

// Default implementation.  Note, shell matrices have to do something else!
void MatrixBase::display() const
{
  MatView(d_A, PETSC_VIEWER_STDOUT_SELF);
}

} // end namespace linear_algebra

//---------------------------------------------------------------------------//
//              end of file MatrixBase.cc
//---------------------------------------------------------------------------//
