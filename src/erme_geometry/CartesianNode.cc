//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   CartesianNode.cc
 *  @brief  CartesianNode
 *  @author Jeremy Roberts
 *  @date   Aug 22, 2012
 */
//---------------------------------------------------------------------------//

#include "CartesianNode.hh"

namespace erme_geometry
{

CartesianNode::CartesianNode(const size_t  d,
                             const size_t  n,
                             const int     nodeid,
                             std::string   nodename,
                             const Point   nodeorigin,
                             vec2_size_t   so,
                             vec_size_t    po,
                             vec_size_t    ao,
                             vec_size_t    eo,
                             vec_dbl          nodewidth)
  : Node(d, n, nodeid, nodename, nodeorigin, so, po, ao, eo)
  , d_width(nodewidth)
{
  // Preconditions
  Require(nodewidth.size() == 3);
  Require(number_surfaces() == dimension()*2);

  // Make sure slabs and boxes get dummy areas and volumes
  if (dimension() < 3) d_width[2] = 1.0;
  if (dimension() < 2) d_width[1] = 1.0;
}

} // end namespace erme_geometry

BOOST_CLASS_EXPORT_IMPLEMENT(erme_geometry::CartesianNode)

//---------------------------------------------------------------------------//
//              end of file CartesianNode.cc
//---------------------------------------------------------------------------//
