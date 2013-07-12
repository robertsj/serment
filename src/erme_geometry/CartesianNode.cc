//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  CartesianNode.cc
 *  @brief CartesianNode member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "CartesianNode.hh"
#include "NodeSerialization.hh"
#include "geometry/Point.hh"
#include <cmath>

namespace erme_geometry
{

//----------------------------------------------------------------------------//
CartesianNode::CartesianNode(const size_t  dim,
                             std::string   nodename,
                             vec2_size_t   so,
                             vec_size_t    po,
                             vec_size_t    ao,
                             vec_size_t    eo,
                             vec_dbl       nodewidth)
  : Node(dim, 2*dim, nodename, so, po, ao, eo)
  , d_width(nodewidth)
{
  Require(nodewidth.size() == 3);
  Require(number_surfaces() == dim*2);

  // Make sure slabs and boxes get dummy areas and volumes
  if (dimension() < 3) d_width[2] = 1.0;
  if (dimension() < 2) d_width[1] = 1.0;
}

//----------------------------------------------------------------------------//
CartesianNode::SP_node
CartesianNode::Create(const size_t  dimension,
                      std::string   nodename,
                      vec2_size_t   so,
                      vec_size_t    po,
                      vec_size_t    ao,
                      vec_size_t    eo,
                      vec_dbl       nodewidth)
{
  SP_node p(new CartesianNode(dimension, nodename,
                              so, po, ao, eo, nodewidth));
  return p;
}

//----------------------------------------------------------------------------//
double CartesianNode::area(const size_t surface) const
{
  Require(surface < number_surfaces());
  if (surface == 0 or surface == 1)
    return d_width[1] * d_width[2];
  if (surface == 2 or surface == 3)
    return d_width[0] * d_width[2];
  return d_width[0] * d_width[1];
}

//----------------------------------------------------------------------------//
double CartesianNode::volume() const
{
  return d_width[0] * d_width[1] * d_width[2];
}

//----------------------------------------------------------------------------//
double CartesianNode::color(const Point &point, std::string)
{
  Point point_local = point;
  if ( (point_local.x() < 0.0) or
       (point_local.y() < 0.0) or
       (point_local.z() < 0.0) or
       (point_local.x() > d_width[0]) or
       (point_local.y() > d_width[1]) or
       (point_local.z() > d_width[2]) )
  {
    return -1.0;
  }
  return (double) 1.0;
}

//----------------------------------------------------------------------------//
double CartesianNode::width(const size_t dim) const
{
  Require(dim < 3);
  return d_width[dim];
}

} // end namespace erme_geometry

BOOST_CLASS_EXPORT_IMPLEMENT(erme_geometry::CartesianNode)

//----------------------------------------------------------------------------//
//              end of file CartesianNode.cc
//----------------------------------------------------------------------------//
