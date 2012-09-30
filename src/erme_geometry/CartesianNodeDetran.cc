//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CartesianNodeDetran.cc
 * \brief  CartesianNodeDetran member definitions
 * \author Jeremy Roberts
 * \date   Aug 22, 2012
 */
//---------------------------------------------------------------------------//

#include "CartesianNodeDetran.hh"

#include "utilities/detran_utilities.hh"

namespace erme_geometry
{

CartesianNodeDetran::CartesianNodeDetran(const size_type  d,
                                         const size_type  n,
                                         const int        nodeid,
                                         std::string      nodename,
                                         const Point      nodeorigin,
                                         vec2_size_type   so,
                                         vec_size_type    po,
                                         vec_size_type    ao,
                                         vec_size_type    eo,
                                         vec_dbl          nodewidth,
                                         SP_db            nodedb,
                                         SP_material      nodematerial,
                                         SP_mesh          nodemesh)
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

#ifdef SERMENT_ENABLE_BOOST
BOOST_CLASS_EXPORT_IMPLEMENT(erme_geometry::CartesianNodeDetran)
#endif

//---------------------------------------------------------------------------//
//              end of file CartesianNodeDetran.cc
//---------------------------------------------------------------------------//
