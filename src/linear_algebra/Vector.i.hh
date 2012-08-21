//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Vector.i.hh
 * \brief  Vector.i 
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 */
//---------------------------------------------------------------------------//

#ifndef VECTOR_I_HH_
#define VECTOR_I_HH_

namespace linear_algebra
{

// Value Setting

inline void Vector::insert_values(const size_type number,
                                  const int *rows,
                                  const double *values)
{
  VecSetValues(d_V, number, rows, values, INSERT_VALUES);
  d_is_assembled = false;
}

//---------------------------------------------------------------------------//
// Vector Operations
//---------------------------------------------------------------------------//

inline double Vector::dot(Vector &x)
{
  // Preconditions
  Require(x.is_assembled());
  Require(is_assembled());
  Require(x.global_size() == d_global_size);

  PetscErrorCode ierr;
  double val;
  ierr = VecDot(d_V, x.V(), &val);

  // Postconditions
  Ensure(!ierr);
  return val;
}

inline void Vector::scale(const double factor)
{
  // Preconditions
  Require(is_assembled());

  PetscErrorCode ierr;
  ierr = VecScale(d_V, factor);

  // Postconditions
  Ensure(!ierr);
}

//---------------------------------------------------------------------------//
// Accessors
//---------------------------------------------------------------------------//

inline const double&
Vector::operator[](const size_type i) const
{
  Require(i < d_local_size);
  return d_array[i];
}

inline double&
Vector::operator[](const size_type i)
{
  // Cast away return type
  return const_cast<double&>
  (
    // Add const to *this's type and call const version
    static_cast<const Vector&>(*this)[i]
  );
}

} // end namespace linear_algebra

#endif // VECTOR_I_HH_ 

//---------------------------------------------------------------------------//
//              end of file Vector.i.hh
//---------------------------------------------------------------------------//
