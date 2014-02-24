//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  NodeResponse.i.hh
 *  @brief NodeResponse inline member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef erme_response_NODERESPONSE_I_HH_
#define erme_response_NODERESPONSE_I_HH_

#include <iostream>

namespace erme_response
{

//----------------------------------------------------------------------------//
inline double
NodeResponse::boundary_response(const size_t out,
                                const size_t in) const
{
  Require(in  < d_N);
  Require(out < d_N);
  return d_boundary_response[in][out];
}
//----------------------------------------------------------------------------//
inline double&
NodeResponse::boundary_response(const size_t out,
                                const size_t in)
{
  Require(in  < d_N);
  Require(out < d_N);
  return d_boundary_response[in][out];
}

//----------------------------------------------------------------------------//
inline double
NodeResponse::fission_response(const size_t in) const
{
  Require(in < d_N);
  return d_fission_response[in];
}
//----------------------------------------------------------------------------//
inline double&
NodeResponse::fission_response(const size_t in)
{
  Require(in < d_N);
  return d_fission_response[in];
}

//----------------------------------------------------------------------------//
inline double
NodeResponse::absorption_response(const size_t in) const
{
  Require(in < d_N);
  return d_absorption_response[in];
}
//----------------------------------------------------------------------------//
inline double&
NodeResponse::absorption_response(const size_t in)
{
  Require(in < d_N);
  return d_absorption_response[in];
}

//----------------------------------------------------------------------------//
inline double
NodeResponse::leakage_response(const size_t surface,
                               const size_t in) const
{
  Require(surface < d_number_surfaces);
  Require(in < d_N);
  return d_leakage_response[in][surface];
}
//----------------------------------------------------------------------------//
inline double&
NodeResponse::leakage_response(const size_t surface,
                               const size_t in)
{
  Require(surface < d_number_surfaces);
  Require(in < d_N);
  return d_leakage_response[in][surface];
}

//----------------------------------------------------------------------------//
inline double
NodeResponse::nodal_power(const size_t in) const
{
  Require(in < d_N);
  return d_nodal_power[in];
}
//----------------------------------------------------------------------------//
inline double&
NodeResponse::nodal_power(const size_t in)
{
  Require(in < d_N);
  return d_nodal_power[in];
}

//----------------------------------------------------------------------------//
inline double
NodeResponse::pin_power(const size_t p,
                        const size_t in) const
{
  Requirev(p < d_number_pins, "p = " + AsString(p));
  Requirev(in < d_N, "in = " + AsString(in));
  return d_pin_power[in][p];
}
//----------------------------------------------------------------------------//
inline double&
NodeResponse::pin_power(const size_t p,
                        const size_t in)
{
  Requirev(p < d_number_pins, "p = " + AsString(p));
  Requirev(in < d_N, "in = " + AsString(in));
  return d_pin_power[in][p];
}

//----------------------------------------------------------------------------//
inline void NodeResponse::clear()
{
	for (size_t i = 0; i < d_N; ++i)
	{
		for (size_t j = 0; j < d_N; ++j)
			d_boundary_response[i][j] = 0.0;
		for (size_t j = 0; j < d_number_surfaces; ++j)
			d_leakage_response[i][j] = 0.0;
		d_fission_response[i] = 0.0;
		d_absorption_response[i] = 0.0;
		d_nodal_power[i] = 0.0;
    for (size_t j = 0; j < d_number_pins; ++j)
      d_pin_power[i][j] = 0.0;
	}
}

} // end namespace erme_response

#endif // erme_response_NODERESPONSE_I_HH_

//----------------------------------------------------------------------------//
//              end of file NodeResponse.i.hh
//----------------------------------------------------------------------------//
