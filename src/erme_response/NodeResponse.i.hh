//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   NodeResponse.i.hh
 * \brief  NodeResponse inline member definitions
 * \author Jeremy Roberts
 * \date   Aug 28, 2012
 */
//---------------------------------------------------------------------------//

#ifndef NODERESPONSE_I_HH_
#define NODERESPONSE_I_HH_

#include <iostream>

namespace erme_response
{

//-------------------------------------------------------------------------//
// ACCESS
//-------------------------------------------------------------------------//

/// Const access to boundary response
inline const double&
NodeResponse::boundary_response(const size_t out,
                                const size_t in) const
{
  Require(in  < d_N);
  Require(out < d_N);
  return d_boundary_response[in][out];
}

/// Mutable access to boundary response
inline double&
NodeResponse::boundary_response(const size_t out,
                                const size_t in)
{
  // Cast away return type
  return const_cast<double&>
  (
    // Add const to *this's type and call const version
    static_cast<const NodeResponse*>(this)->boundary_response(out, in)
  );
}

/// Const access to fission response
inline const double&
NodeResponse::fission_response(const size_t in) const
{
  Require(in < d_N);
  return d_fission_response[in];
}

/// Mutable access to fission response
inline double&
NodeResponse::fission_response(const size_t in)
{
  // Cast away return type
  return const_cast<double&>
  (
    // Add const to *this's type and call const version
    static_cast<const NodeResponse*>(this)->fission_response(in)
  );
}

/// Const access to absorption response
inline const double&
NodeResponse::absorption_response(const size_t in) const
{
  Require(in < d_N);
  return d_absorption_response[in];
}

/// Mutable access to absorption response
inline double&
NodeResponse::absorption_response(const size_t in)
{
  // Cast away return type
  return const_cast<double&>
  (
    // Add const to *this's type and call const version
    static_cast<const NodeResponse*>(this)->absorption_response(in)
  );
}

/// Const access to leakage response
inline const double&
NodeResponse::leakage_response(const size_t surface,
                               const size_t in) const
{
  Require(surface < d_number_surfaces);
  Require(in < d_N);
  return d_leakage_response[in][surface];
}

/// Mutable access to leakage response
inline double&
NodeResponse::leakage_response(const size_t surface,
                               const size_t in)
{
  // Cast away return type
  return const_cast<double&>
  (
    // Add const to *this's type and call const version
    static_cast<const NodeResponse*>(this)->leakage_response(surface, in)
  );
}

} // end namespace erme_response

#endif // NODERESPONSE_I_HH_ 

//---------------------------------------------------------------------------//
//              end of file NodeResponse.i.hh
//---------------------------------------------------------------------------//
