//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Vector.hh
 *  @brief Vector class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef linear_algebra_VECTOR_HH_
#define linear_algebra_VECTOR_HH_

#include "utilities/DBC.hh"
#include "utilities/SP.hh"
#include "petsc.h"

namespace linear_algebra
{

/**
 *  @class Vector
 *  @brief Lightweight wrapper for PETSc Vec.
 *
 *  Use of Vector and the corresponding \ref Matrix class should, in
 *  theory, eliminate a lot of PETSc code from the rest of Serment.
 *
 */
/**
 *  @example linear_algebra/test/test_Vector.cc
 *
 *  Test of Vector class.
 */
class Vector
{

public:

  //--------------------------------------------------------------------------//
  // ENUMERATIONS
  //--------------------------------------------------------------------------//

  /**
   *  @brief various vector norms available
   *
   *  For a vector @f$ v \in (m, 1) @f$:
   *    - @f$ ||v||_1         \equiv \sum_i |v_i| @f$
   *    - @f$ ||v||_2         \equiv \sqrt{\sum_i v_i^2} @f$
   *    - @f$ ||v||_{\infty}  \equiv \max_{v_i} |v_i| @f$
   */
  enum vector_norm_types
  {
    L1, L2, LINF, END_VECTOR_NORM_TYPES
  };


  enum vector_display_types
  {
    STDOUT, ASCII, BINARY, END_VECTOR_DISPLAY_TYPES
  };

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<Vector>      SP_vector;
  typedef detran_utilities::size_t          size_type;
  typedef detran_utilities::vec_int         vec_int;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param m      Local number of rows
   *  @param val    Optional initial value
   */
  Vector(const size_type m, const double val = 0.0, const bool seq = false);

  /**
   *  @brief Copy Constructor
   *  @param V      Vector to copy
   */
  Vector(const Vector &V);

  /// Temporarily wrap a vector with a smaller vector: does *NOT* free memory
  Vector(Vector &V, const int m);

  /// Temporary construction via PETSc vector: does *NOT* free memory
  Vector(Vec pv);

  /// Destructor
  ~Vector();

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /**
   *  @brief Insert values
   *  @param values   Array of values to insert
   *  @param number   Number of values to insert
   *  @param rows     Indices of rows where values are inserted
   */
  void insert_values(const unsigned int number,
                     const int *rows,
                     const double *values);


  /// Assemble the vector.
  void assemble();

  //@{
  /// Vector operations
  double norm(const int type = L2);
  double norm_residual(const Vector& x, const int type = L2);
  double dot(Vector &x);
  void scale(const double factor);
  void set(const double v);
  void add(const Vector& x);
  void subtract(const Vector& x);
  void multiply(const Vector& x);
  void divide(const Vector& x);
  void copy(const Vector& x);
  void add_a_times_x(const double a, const Vector& x);
  //@}

  ///@{
  /// Access to element
  const double& operator[](const size_type i) const;
  double& operator[](const size_type i);
  ///@}

  //@{
  /// Getters
  Vec V() const;
  int global_size() const;
  int local_size() const;
  int lower_bound() const;
  int upper_bound() const;
  vec_int ranges() const;
  bool is_assembled() const;
  bool is_temporary() const;
  bool is_sequential() const;
  //@}

  /**
   *  @brief Collect a multi-process vector on a single vector
   *  @param    v     Parallel vector
   *  @param    root  Process on which sequential vector is stored
   *  @return         Sequential vector containing parallel vector contents
   */
  SP_vector collect_on_root(const size_type root = 0);


  /// View to screen, to ascii file, or to binary file
  void display(const size_type    display_type = STDOUT,
               const std::string &filename = "V.out") const;

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

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
  /// Am I temporary?
  bool d_is_temporary;
  /// Am I a sequential vector?
  bool d_is_sequential;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  /// Consolidation of common construction steps
  PetscErrorCode setup();

};

} // end namespace linear_algebra

//----------------------------------------------------------------------------//
// INLINE MEMBER DEFINITIONS
//----------------------------------------------------------------------------//

#include "Vector.i.hh"

#endif // linear_algebra_VECTOR_HH_

//----------------------------------------------------------------------------//
//              end of file Vector.hh
//----------------------------------------------------------------------------//
