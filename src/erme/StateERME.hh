//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   StateERME.hh
 * \brief  StateERME 
 * \author Jeremy Roberts
 * \date   Aug 23, 2012
 */
//---------------------------------------------------------------------------//

#ifndef STATEERME_HH_
#define STATEERME_HH_

#include "linear_algebra/Vector.hh"

namespace erme
{

class StateERME
{

public:

  typedef linear_algebra::Vector moments_type;
  typedef unsigned int           size_type;

  /*!
   *  \brief Constructor
   *  \param size   Size of local state vector
   */
  StateERME(const int size_type);

  // SETTERS

  void set_k(const double k_val);

  void set_lambda(const double lambda_val);

  // GETTERS

  double k() const;

  double lambda() const;

  // MOMENT ACCESS



private:

  /// Boundary unknowns
  moments_type d_boundary_moments;

  /// K-eigenvalue
  double d_k;

  /// Lambda-eigenvalue
  double d_lambda;

};

} // end namespace erme

// Inline member definitions
#include "StateERME.i.hh"

#endif // STATEERME_HH_ 

//---------------------------------------------------------------------------//
//              end of file StateERME.hh
//---------------------------------------------------------------------------//
