//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Vector.hh
 * \brief  Vector 
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 */
//---------------------------------------------------------------------------//

#ifndef VECTOR_HH_
#define VECTOR_HH_

// Detran
#include "DBC.hh"

// System
#include "petsc.h"

namespace linear_algebra
{

/*!
 *  \class Vector
 *  \brief Lightweight wrapper for PETSc Vec.
 *
 *  Use of Vector and the corresponding \ref Matrix class should, in
 *  theory, eliminate a lot of PETSc code from the rest of Serment.
 *
 */
/*!
 *  \example linear_algebra/test/test_Vector.cc
 *
 *  Test of Vector class.
 */
class Vector
{

public:

  typedef unsigned int size_type;

  /*!
   *  \brief Constructor
   *  \param m      Local number of rows
   *  \param val    Optional initial value
   */
  Vector(const size_type m, const double val = 0.0);

  /// Destructor
  ~Vector();

  //---------------------------------------------------------------------------//
  // SETTERS
  //---------------------------------------------------------------------------//

  /*!
   *  \brief Insert values
   *  \param values   Array of values to insert
   *  \param number   Number of values to insert
   *  \param rows     Indices of rows where values are inserted
   */
  void insert_values(const unsigned int number,
                     const int *rows,
                     const double *values);


  /// Assemble the vector.
  void assemble();

  //---------------------------------------------------------------------------//
  // VECTOR OPERATIONS
  //---------------------------------------------------------------------------//

  /// Dot product of another Vector with me
  double dot(Vector &x);

  /// Scale the Vector
  void scale(const double factor);

  //---------------------------------------------------------------------------//
  // ACCESSORS
  //---------------------------------------------------------------------------//

  /*
   *  \brief Const access to local array
   *  \param i  Local index
   */
  const double& operator[](const size_type i) const;

  /*
   *  \brief Mutable access to local array
   *  \param i  Local index
   */
  double& operator[](const size_type i);

  /// Return the PETSc vector
  Vec V() const
  {
    return d_V;
  }

  //---------------------------------------------------------------------------//
  // GETTERS
  //---------------------------------------------------------------------------//

  /// Return the global size.
  int global_size() const
  {
    return d_global_size;
  }

  /// Return the local size.
  int local_size() const
  {
    return d_local_size;
  }

  /// Return the lower bound.
  int lower_bound() const
  {
    return d_lower_bound;
  }

  /// Return the upper bound.
  int upper_bound() const
  {
    return d_upper_bound;
  }

  /// Return assembled flag.
  bool is_assembled() const
  {
    return d_is_assembled;
  }

  /// View via standard output.
  void display() const;

private:

  //---------------------------------------------------------------------------//
  // PRIVATE DATA
  //---------------------------------------------------------------------------//

  /// PETSc vector
  Vec d_V;

  /// Pointer to underlying local array.  Be *careful* when using.
  double* d_array;

  /// Global size
  size_type d_global_size;

  /// Local size
  size_type d_local_size;

  /// Local lower bound
  size_type d_lower_bound;

  /// Local upper bound
  size_type d_upper_bound;

  /// Am I assembled?
  bool d_is_assembled;

};

} // end namespace linear_algebra

// Inline members
#include "Vector.i.hh"

#endif // VECTOR_HH_ 

//---------------------------------------------------------------------------//
//              end of file Vector.hh
//---------------------------------------------------------------------------//
