//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  StateERME.cc
 *  @brief StateERME member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "StateERME.hh"

namespace erme
{

//----------------------------------------------------------------------------//
StateERME::StateERME(const size_t size)
  : d_local_size(0)
  , d_k(1.0)
  , d_lambda(1.0)
{
  if (serment_comm::Comm::is_global())
  {
  	d_boundary_moments = new Vector(size, 0.0);
    d_global_size = d_boundary_moments->global_size();
    d_local_size = size;
  }
  serment_comm::Comm::broadcast(&d_global_size, 1, 0);
}

//----------------------------------------------------------------------------//
void StateERME::update(SP_vector v, const double k, const double l)
{
	if (serment_comm::Comm::is_global())
	{
		Require(v);
		(*d_boundary_moments).copy(*v);
	}
  d_k = k;
  d_lambda = l;
  serment_comm::Comm::broadcast(&d_k, 1, 0);
  serment_comm::Comm::broadcast(&d_lambda, 1, 0);
}

} // end namespace erme

//----------------------------------------------------------------------------//
//              end of file StateERME.cc
//----------------------------------------------------------------------------//
