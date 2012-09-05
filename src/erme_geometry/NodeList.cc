//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   NodeList.cc
 * \brief  NodeList member definitions
 * \author Jeremy Roberts
 * \date   Aug 27, 2012
 */
//---------------------------------------------------------------------------//

#include "NodeList.hh"

namespace erme_geometry
{

NodeList::NodeList()
  : d_lower_bound(0),
    d_upper_bound(0),
    d_number_local_surfaces(0),
    d_number_global_surfaces(0),
    d_is_finalized(false)
{}

} // end namespace erme_geometry

//---------------------------------------------------------------------------//
//              end of file NodeList.cc
//---------------------------------------------------------------------------//
