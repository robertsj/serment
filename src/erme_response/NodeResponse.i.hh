//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   NodeResponse.i.hh
 *  @brief  NodeResponse inline member definitions
 *  @author Jeremy Roberts
 *  @date   Aug 28, 2012
 */
//---------------------------------------------------------------------------//

#ifndef erme_response_NODERESPONSE_I_HH_
#define erme_response_NODERESPONSE_I_HH_

#include <iostream>

namespace erme_response
{

//-------------------------------------------------------------------------//
inline const double&
NodeResponse::boundary_response(const size_t out,
                                const size_t in) const
{
  Require(in  < d_N);
  Require(out < d_N);
  return d_boundary_response[in][out];
}
//-------------------------------------------------------------------------//
inline double&
NodeResponse::boundary_response(const size_t out,
                                const size_t in)
{
  return const_cast<double&>
  (
    static_cast<const NodeResponse*>(this)->boundary_response(out, in)
  );
}

//-------------------------------------------------------------------------//
inline const double&
NodeResponse::fission_response(const size_t in) const
{
  Require(in < d_N);
  return d_fission_response[in];
}
//-------------------------------------------------------------------------//
inline double&
NodeResponse::fission_response(const size_t in)
{
  return const_cast<double&>
  (
    static_cast<const NodeResponse*>(this)->fission_response(in)
  );
}

//-------------------------------------------------------------------------//
inline const double&
NodeResponse::absorption_response(const size_t in) const
{
  Require(in < d_N);
  return d_absorption_response[in];
}
//-------------------------------------------------------------------------//
inline double&
NodeResponse::absorption_response(const size_t in)
{
  return const_cast<double&>
  (
    static_cast<const NodeResponse*>(this)->absorption_response(in)
  );
}

//-------------------------------------------------------------------------//
inline const double&
NodeResponse::leakage_response(const size_t surface,
                               const size_t in) const
{
  Require(surface < d_number_surfaces);
  Require(in < d_N);
  return d_leakage_response[in][surface];
}
//-------------------------------------------------------------------------//
inline double&
NodeResponse::leakage_response(const size_t surface,
                               const size_t in)
{
  return const_cast<double&>
  (
    static_cast<const NodeResponse*>(this)->leakage_response(surface, in)
  );
}

} // end namespace erme_response

#endif // erme_response_NODERESPONSE_I_HH_

//---------------------------------------------------------------------------//
//              end of file NodeResponse.i.hh
//---------------------------------------------------------------------------//
