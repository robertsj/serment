//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   StateERME.cc
 * \brief  StateERME member definitions
 * \author Jeremy Roberts
 * \date   Aug 23, 2012
 */
//---------------------------------------------------------------------------//

#include "StateERME.hh"

namespace erme
{

StateERME::StateERME(const int size)
  : d_boundary_moments(size)
  , d_k(1.0)
  , d_lambda(1.0)
{
  /* ... */
}

} // end namespace erme

//---------------------------------------------------------------------------//
//              end of file StateERME.cc
//---------------------------------------------------------------------------//
