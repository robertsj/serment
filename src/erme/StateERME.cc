//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   StateERME.cc
 *  @brief  StateERME member definitions
 *  @author Jeremy Roberts
 *  @date   Aug 23, 2012
 */
//---------------------------------------------------------------------------//

#include "StateERME.hh"

namespace erme
{

//---------------------------------------------------------------------------//
StateERME::StateERME(const size_t size)
  : d_boundary_moments(size, 0.0)
  , d_local_size(size)
  , d_k(1.0)
  , d_lambda(1.0)
{
  Require(size == d_boundary_moments.local_size());
  d_global_size = d_boundary_moments.global_size();
}

} // end namespace erme

//---------------------------------------------------------------------------//
//              end of file StateERME.cc
//---------------------------------------------------------------------------//
