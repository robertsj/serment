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
 *  \brief Lightweight wrapper for PETSc Vec
 *
 */
class Vector
{

public:

  /*!
   *  \brief Constructor
   *  \param m    Number of rows
   */
  Vector(const unsigned int m);

  /// Destructor
  ~Vector();


  /// \name Setters
  /// \{

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

  /// \}

  /// \name Vector Operators
  /// \{

  /// Dot product of another Vector with me
  double dot(Vector &x);

  /// Scale the Vector
  void scale(const double factor);

  /// \}

  /// Return the PETSc vector
  Vec V()
  {
    return d_V;
  }

  /// Return the global size
  int global_size() const
  {
    return d_global_size;
  }

  /// Return the local size
  int local_size() const
  {
    return d_local_size;
  }

  /// View via standard or binary output
  void display() const;

  /// Return assembled flag
  bool is_assembled() const
  {
    return d_is_assembled;
  }

private:

  /// \name Private Data
  /// \{

  /// PETSc vector
  Vec d_V;

  /// Global size
  unsigned int d_global_size;

  /// Local size
  unsigned int d_local_size;

  /// Am I assembled?
  bool d_is_assembled;

  /// \}

};


} // end namespace linear_algebra

// Inline members
#include "Vector.i.hh"

#endif // VECTOR_HH_ 

//---------------------------------------------------------------------------//
//              end of file Vector.hh
//---------------------------------------------------------------------------//
