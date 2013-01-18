//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   CartesianNodeDetran.cc
 *  @brief  CartesianNodeDetran member definitions
 *  @author Jeremy Roberts
 *  @date   Aug 22, 2012
 */
//---------------------------------------------------------------------------//

#include "CartesianNodeDetran.hh"
#include "NodeSerialization.hh"

#include "utilities/detran_utilities.hh"

namespace erme_geometry
{

//---------------------------------------------------------------------------//
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
  Require(d_db);
  Require(d_material);
  Require(d_mesh);
  Require(d_mesh->dimension() == dimension());
  Require(detran_utilities::soft_equiv(d_mesh->total_width_x(), width(0)));
  Require(detran_utilities::soft_equiv(d_mesh->total_width_y(), width(1)));
  Require(detran_utilities::soft_equiv(d_mesh->total_width_z(), width(2)));
}

//---------------------------------------------------------------------------//
double CartesianNodeDetran::color(Point point, std::string key)
{
  size_t ijk[] = {0, 0, 0};
  double xyz[] = {point.x(), point.y(), point.z()};
  for (size_t d = 0; d < d_mesh->dimension(); ++d)
  {
    double uvw = 0.0;
    size_t i;
    for (i = 0; i < d_mesh->number_cells(d); ++i)
      if (xyz[d] < uvw + d_mesh->width(d, i)) break;
    ijk[d] = i;
  }
  size_t idx = d_mesh->index(ijk[0], ijk[1], ijk[2]);
  Insist(d_mesh->mesh_map_exists(key),
         "Key used for CartesianNodeDetran coloring does not exist!");
  const vec_int &mat_map = d_mesh->mesh_map(key);
  return (double)mat_map[idx];
}

} // end namespace erme_geometry

BOOST_CLASS_EXPORT_IMPLEMENT(erme_geometry::CartesianNodeDetran)

//---------------------------------------------------------------------------//
//              end of file CartesianNodeDetran.cc
//---------------------------------------------------------------------------//
