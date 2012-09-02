//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   NodeResponse.cc
 * \brief  NodeResponse member definitions
 * \author Jeremy Roberts
 * \date   Aug 28, 2012
 */
//---------------------------------------------------------------------------//

#include "NodeResponse.hh"

namespace erme_response
{

NodeResponse::NodeResponse(const size_t N, const size_t number_surfaces)
  : d_N(N)
  , d_number_surfaces(number_surfaces)
  , d_boundary_response(N, vec_dbl(N, 0.0))
  , d_fission_response(N, 0.0)
  , d_absorption_response(N, 0.0)
  , d_leakage_response(N, vec_dbl(number_surfaces, 0.0))
{
  Require(d_N > 1); // At least 2 surfaces, each with at least 1 moment
  Require(d_number_surfaces > 1);
}

} // end namespace erme_response

//---------------------------------------------------------------------------//
//              end of file NodeResponse.cc
//---------------------------------------------------------------------------//
