//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Node.cc
 * \brief  Node 
 * \author Jeremy Roberts
 * \date   Aug 22, 2012
 */
//---------------------------------------------------------------------------//

#include "Node.hh"

namespace erme_geometry
{

Node::Node(const size_type  d,
           const size_type  n,
           const int        nodeid,
           std::string      nodename,
           const Point      nodeorigin,
           vec2_size_type   so,
           vec_size_type    po,
           vec_size_type    ao,
           vec_size_type    eo)
  : d_dimension(d)
  , d_number_surfaces(n)
  , d_id(nodeid)
  , d_name(nodename)
  , d_origin(nodeorigin)
  , d_spatial_order(so)
  , d_polar_order(po)
  , d_azimuthal_order(ao)
  , d_energy_order(eo)
{
  // Preconditions
  Require(d_dimension > 0 and d_dimension <= 3);
  Require(d_number_surfaces > 0);
  Require(d_spatial_order.size()    == d_number_surfaces);
  Require(d_polar_order.size()      == d_number_surfaces);
  Require(d_azimuthal_order.size()  == d_number_surfaces);
  Require(d_energy_order.size()     == d_number_surfaces);
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file Node.cc
//---------------------------------------------------------------------------//