//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   StateERME.hh
 *  @brief  StateERME
 *  @author Jeremy Roberts
 *  @date   Aug 23, 2012
 */
//---------------------------------------------------------------------------//

#ifndef erme_STATEERME_HH_
#define erme_STATEERME_HH_

#include "linear_algebra/Vector.hh"

namespace erme
{

/**
 *  @class State
 *  @brief Represents the problem state vector
 *
 *  A solution for the eigenvalue response matrix equations consists
 *  of a global boundary vector, which contains moments for each surface
 *  of each cell, and the k-eigenvalue.
 *
 */
/**
 *  @example erme/test/test_StateERME.cc
 *
 *  Test of StateERME class
 */
class StateERME
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<StateERME>   SP_state;
  typedef linear_algebra::Vector            Vector;
  typedef unsigned int                      size_t;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param size   Size of local state vector
   */
  StateERME(const size_t size);

  // SETTERS
  void update(const Vector &v, const double k, const double l);
  void set_k(const double k_val);
  void set_lambda(const double lambda_val);

  // GETTERS
  double k() const;
  double lambda() const;
  size_t local_size() const;
  size_t global_size() const;

  // MOMENT ACCESS (for now, just return the vector directly)

  /// Const reference to moments vector
  const Vector& moments() const;
  /// Mutable reference to moments vector
  Vector& moments();

private:

  /// Boundary unknowns
  Vector d_boundary_moments;
  /// Local size of moments
  size_t d_local_size;
  /// Global size of moments
  size_t d_global_size;
  /// K-eigenvalue
  double d_k;
  /// Lambda-eigenvalue
  double d_lambda;

};

} // end namespace erme

// Inline member definitions
#include "StateERME.i.hh"

#endif // erme_STATEERME_HH_

//---------------------------------------------------------------------------//
//              end of file StateERME.hh
//---------------------------------------------------------------------------//
