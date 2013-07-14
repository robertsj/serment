//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Vector.i.hh
 *  @brief Vector inline member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef linear_algebra_VECTOR_I_HH_
#define linear_algebra_VECTOR_I_HH_

namespace linear_algebra
{

//----------------------------------------------------------------------------//
inline void Vector::insert_values(const size_type  number,
                                  const int       *rows,
                                  const double    *values)
{
  VecSetValues(d_V, number, rows, values, INSERT_VALUES);
  d_is_assembled = false;
}

//----------------------------------------------------------------------------//
inline double Vector::norm(const int type)
{
  Require(is_assembled());

  double val = 0.0;
  PetscErrorCode ierr;
  if (type == L2)
    ierr = VecNorm(d_V, NORM_2, &val);
  else if (type == L1)
    ierr = VecNorm(d_V, NORM_1, &val);
  else if (type == LINF)
    ierr = VecNorm(d_V, NORM_INFINITY, &val);
  else
    THROW("Unsupported vector norm type");
  Ensure(!ierr);
  Ensure(val >= 0.0);//, "Invalid norm: " + AsString(val));
  return val;
}

//----------------------------------------------------------------------------//
inline double Vector::norm_residual(const Vector &x, const int type)
{
  Require(x.is_assembled());
  Require(is_assembled());
  Require(x.global_size() == d_global_size);

  // Create residual vector
  Vector r(x.local_size(), 0.0);
  // Fill it
  r.add((*this));
  r.subtract(x);

  return r.norm(type);
}

//----------------------------------------------------------------------------//
inline double Vector::dot(Vector &x)
{
  Require(x.is_assembled());
  Require(is_assembled());
  Require(x.global_size() == d_global_size);

  PetscErrorCode ierr;
  double val;
  ierr = VecDot(d_V, x.V(), &val);

  if (ierr)
  {
    std::cout << " BAD VECTORS!!! " << std::endl;
    this->display();
    x.display();
  }

  Ensure(!ierr);
  return val;
}

//----------------------------------------------------------------------------//
inline void Vector::scale(const double factor)
{
  Require(is_assembled());

  PetscErrorCode ierr;
  ierr = VecScale(d_V, factor);

  Ensure(!ierr);
}

//----------------------------------------------------------------------------//
inline void Vector::set(const double factor)
{
  Require(is_assembled());

  PetscErrorCode ierr;
  ierr = VecSet(d_V, factor);

  Ensure(!ierr);
}

//----------------------------------------------------------------------------//
inline void Vector::add(const Vector& x)
{
  Require(is_assembled());

  PetscErrorCode ierr;
  ierr = VecAXPY(d_V, 1.0, x.V());

  Ensure(!ierr);
}

//----------------------------------------------------------------------------//
inline void Vector::subtract(const Vector& x)
{
  Require(is_assembled());

  PetscErrorCode ierr;
  ierr = VecAXPY(d_V, -1.0, x.V());

  Ensure(!ierr);
}

//----------------------------------------------------------------------------//
inline void Vector::add_a_times_x(const double a, const Vector& x)
{
  Require(x.local_size() == local_size());

  PetscErrorCode ierr;
  ierr = VecAXPY(d_V, a, const_cast<Vector*>(&x)->V());

  Ensure(!ierr);
}

//----------------------------------------------------------------------------//
inline void Vector::multiply(const Vector& x)
{
  Require(is_assembled());

  // Copy myself
  Vector tmp(*this);

  PetscErrorCode ierr =
    VecPointwiseMult(d_V, tmp.V(), const_cast<Vector*>(&x)->V());

  Ensure(!ierr);
}

//----------------------------------------------------------------------------//
inline void Vector::divide(const Vector& x)
{
  Require(is_assembled());

  // Copy myself
  Vector tmp(*this);

  PetscErrorCode ierr =
    VecPointwiseDivide(d_V, tmp.V(), const_cast<Vector*>(&x)->V());

  Ensure(!ierr);
}

//----------------------------------------------------------------------------//
inline void Vector::copy(const Vector& x)
{
  Require(is_assembled());
  Require(x.is_assembled());
  Require(x.local_size() == local_size());

  PetscErrorCode ierr = VecCopy(x.V(), d_V);

  Ensure(!ierr);
}

//----------------------------------------------------------------------------//
inline const double&
Vector::operator[](const size_type i) const
{
  Require(i < d_local_size);
  return d_array[i];
}

//----------------------------------------------------------------------------//
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


//----------------------------------------------------------------------------//
inline Vec Vector::V() const
{
  return d_V;
}

//----------------------------------------------------------------------------//
inline int Vector::global_size() const
{
  return d_global_size;
}

//----------------------------------------------------------------------------//
inline int Vector::local_size() const
{
  return d_local_size;
}

//----------------------------------------------------------------------------//
inline int Vector::lower_bound() const
{
  return d_lower_bound;
}

//----------------------------------------------------------------------------//
inline int Vector::upper_bound() const
{
  return d_upper_bound;
}

//----------------------------------------------------------------------------//
inline bool Vector::is_assembled() const
{
  return d_is_assembled;
}

//----------------------------------------------------------------------------//
inline bool Vector::is_temporary() const
{
  return d_is_temporary;
}

//----------------------------------------------------------------------------//
inline bool Vector::is_sequential() const
{
  return d_is_sequential;
}

} // end namespace linear_algebra

#endif // linear_algebra_VECTOR_I_HH_

//----------------------------------------------------------------------------//
//              end of file Vector.i.hh
//----------------------------------------------------------------------------//
