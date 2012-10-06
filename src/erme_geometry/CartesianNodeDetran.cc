//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CartesianNodeDetran.cc
 * \brief  CartesianNodeDetran member definitions
 * \author Jeremy Roberts
 * \date   Aug 22, 2012
 */
//---------------------------------------------------------------------------//

#include "CartesianNodeDetran.hh"
#include "NodeSerialization.hh"

#include "utilities/detran_utilities.hh"

namespace erme_geometry
{

CartesianNodeDetran::CartesianNodeDetran(const size_t  d,
                                         const size_t  n,
                                         const int     nodeid,
                                         std::string   nodename,
                                         const Point   nodeorigin,
                                         vec2_size_t   so,
                                         vec_size_t    po,
                                         vec_size_t    ao,
                                         vec_size_t    eo,
                                         vec_dbl       nodewidth,
                                         SP_db         nodedb,
                                         SP_material   nodematerial,
                                         SP_mesh       nodemesh)
  : CartesianNode(d, n, nodeid, nodename, nodeorigin, so, po, ao, eo, nodewidth)
  , d_db(nodedb)
  , d_material(nodematerial)
  , d_mesh(nodemesh)
{
  // Preconditions
  Require(d_db);
  Require(d_material);
  Require(d_mesh);
  Require(d_mesh->dimension() == dimension());
}

// Do something better later.
double CartesianNodeDetran::color(Point point)
{
  return id();
}

} // end namespace erme_geometry

BOOST_CLASS_EXPORT_IMPLEMENT(erme_geometry::CartesianNodeDetran)

//---------------------------------------------------------------------------//
//              end of file CartesianNodeDetran.cc
//---------------------------------------------------------------------------//
