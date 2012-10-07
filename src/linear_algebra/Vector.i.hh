//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Vector.i.hh
 *  @brief  Vector inline member definitions
 *  @author Jeremy Roberts
 *  @date   Aug 19, 2012
 */
//---------------------------------------------------------------------------//

#ifndef VECTOR_I_HH_
#define VECTOR_I_HH_

namespace linear_algebra
{

inline Vector::Vector(const Vector& V)
{
  // Preconditions
  Require(V.is_assembled());

  d_local_size  = V.local_size();
  d_global_size = V.global_size();

  // Create the vector
  PetscErrorCode ierr;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, d_local_size, PETSC_DETERMINE, &d_V);

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

  // Because we initialize to zero, we can assemble right
  // away, since we might just be filling this from a multiplication.
  assemble();

  // Copy the incoming vector
  copy(V);

  // Postconditions
  Ensure(!ierr);
  Ensure(d_upper_bound - d_lower_bound == d_local_size);

}

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

inline double Vector::norm(const int type)
{
  // Preconditions
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

  // Posconditions
  Ensure(!ierr);
  Ensure(val >= 0.0);
  return val;
}

inline double Vector::norm_residual(const Vector &x, const int type)
{
  // Preconditions
  Require(x.is_assembled());
  Require(is_assembled());
  Require(x.global_size() == d_global_size);

  // Create residual vector
  Vector r(x.local_size(), 0.0);
  // Fill it
  r.add((*this));
  r.subtract(x);

  // Return the norm of the residual
  return r.norm(type);
}

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

inline void Vector::set(const double factor)
{
  // Preconditions
  Require(is_assembled());

  PetscErrorCode ierr;
  ierr = VecSet(d_V, factor);

  // Postconditions
  Ensure(!ierr);
}

inline void Vector::add(const Vector& x)
{
  // Preconditions
  Require(is_assembled());

  PetscErrorCode ierr;
  ierr = VecAXPY(d_V, 1.0, x.V());

  // Postconditions
  Ensure(!ierr);
}

inline void Vector::subtract(const Vector& x)
{
  // Preconditions
  Require(is_assembled());

  PetscErrorCode ierr;
  ierr = VecAXPY(d_V, -1.0, x.V());

  // Postconditions
  Ensure(!ierr);
}

inline void Vector::add_a_times_x(const double a, const Vector& x)
{
  // Preconditions
  Require(x.local_size() == local_size());

  PetscErrorCode ierr;
  ierr = VecAXPY(d_V, a, const_cast<Vector*>(&x)->V());

  // Postconditions
  Ensure(!ierr);
}

inline void Vector::multiply(const Vector& x)
{
  // Preconditions
  Require(is_assembled());

  // Copy myself
  Vector tmp(*this);

  PetscErrorCode ierr =
    VecPointwiseMult(d_V, tmp.V(), const_cast<Vector*>(&x)->V());

  // Postconditions
  Ensure(!ierr);
}

inline void Vector::divide(const Vector& x)
{
  // Preconditions
  Require(is_assembled());

  // Copy myself
  Vector tmp(*this);

  PetscErrorCode ierr =
    VecPointwiseDivide(d_V, tmp.V(), const_cast<Vector*>(&x)->V());

  // Postconditions
  Ensure(!ierr);
}

inline void Vector::copy(const Vector& x)
{
  // Preconditions
  Require(is_assembled());
  Require(x.is_assembled());
  Require(x.local_size() == local_size());

  PetscErrorCode ierr = VecCopy(x.V(), d_V);

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
