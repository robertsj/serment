//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  CartesianNodeDetran.cc
 *  @brief CartesianNodeDetran member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "CartesianNodeDetran.hh"
#include "NodeSerialization.hh"
#include "geometry/Point.hh"
#include "utilities/MathUtilities.hh"

#include "utilities/detran_utilities.hh"

namespace erme_geometry
{

//----------------------------------------------------------------------------//
CartesianNodeDetran::CartesianNodeDetran(const size_t  dim,
                                         std::string   nodename,
                                         vec2_size_t   so,
                                         vec_size_t    po,
                                         vec_size_t    ao,
                                         vec_size_t    eo,
                                         vec_dbl       nodewidth,
                                         SP_db         nodedb,
                                         SP_material   nodematerial,
                                         SP_mesh       nodemesh)
  : CartesianNode(dim, nodename, so, po, ao, eo, nodewidth)
  , d_db(nodedb)
  , d_material(nodematerial)
  , d_mesh(nodemesh)
{
  using namespace detran_utilities;

  Require(d_db);
  Require(d_material);
  Require(d_mesh);
  Require(d_mesh->dimension() == dimension());
  Require(soft_equiv(d_mesh->total_width_x(), width(0)));
  Require(soft_equiv(d_mesh->total_width_y(), width(1)));
  Require(soft_equiv(d_mesh->total_width_z(), width(2)));

  if (!d_mesh->mesh_map_exists("NODAL"))
  {
    vec_int tmp(d_mesh->number_cells(), 0);
    d_mesh->add_mesh_map("NODAL", tmp);
  }
  if (d_mesh->mesh_map_exists("PINS"))
  {
    d_number_pins = vec_max(d_mesh->mesh_map("PINS")) + 1;
  }
}

//----------------------------------------------------------------------------//
double CartesianNodeDetran::color(const Point &point, std::string key)
{
  Insist(d_mesh->mesh_map_exists(key),
         "Key used for CartesianNodeDetran coloring does not exist!");

  int cell = d_mesh->find_cell(point);
  if (cell == -1)
  {
    return -1.0;
  }
  else
  {
    const vec_int &mat_map = d_mesh->mesh_map(key);
    return (double)mat_map[cell];
  }
}

} // end namespace erme_geometry

BOOST_CLASS_EXPORT_IMPLEMENT(erme_geometry::CartesianNodeDetran)

//----------------------------------------------------------------------------//
//              end of file CartesianNodeDetran.cc
//----------------------------------------------------------------------------//
