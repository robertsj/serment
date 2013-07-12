//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   node_fixture.hh
 *  @brief  Fixtures for various node types.
 *  @author Jeremy Roberts
 *  @date   Aug 22, 2012
 *
 *  Ideally, fixtures for all available node types should be added to this
 *  file.  Doing so greatly facilitates testing of each node within the
 *  erme_geometry subsystem and Serment as a whole.
 *
 *  Each node should be as simple to compute as possible
 *  (i.e. small, few group, low orders) without limiting the scope of testing.
 *
 */
//----------------------------------------------------------------------------//

#ifndef NODE_FIXTURE_HH_
#define NODE_FIXTURE_HH_

#include "erme_geometry/CartesianNodeDetran.hh"

namespace erme_geometry
{

/**
 *  The Detran Cartesian nodes are all 10 cm in extent
 *  with a monoenergetic material.  All applicable
 *  responses are first order.
 */

Node::SP_node cartesian_node_detran(const int dim,
                                    const int so = 4,
                                    const int ao = 2,
                                    const int po = 2)
{
  Require(dim > 0 and dim <= 3);

  typedef erme_geometry::CartesianNodeDetran Node_T;

  // Parameter database
  Node_T::SP_db   db(new detran_utilities::InputDB());
  db->put("number_groups",  1);
  db->put("dimension",      dim);
  db->put("equation",       "diffusion");

  // Build the material
  Node_T::SP_material mat(new detran_material::Material(1, 1));
  mat->set_sigma_t(0, 0, 1.0);
  mat->set_sigma_s(0, 0, 0, 0.5);
  mat->set_sigma_f(0, 0, 0.5);
  mat->set_sigma_a(0, 0, 0.5);
  mat->finalize();

  // Define discretization and material map
  Node_T::vec_dbl cm(2, 0.0);
  cm[1] = 10.0;
  Node_T::vec_int fm(1, 10);
  Node_T::vec_int mt(1, 0);

  // Build the mesh
  Node_T::SP_mesh mesh;
  if (dim == 1)
    mesh = new detran_geometry::Mesh1D(fm, cm, mt);
  else if (dim == 2)
    mesh = new detran_geometry::Mesh2D(fm, fm, cm, cm, mt);
  else
    mesh = new detran_geometry::Mesh3D(fm, fm, fm, cm, cm, cm, mt);

  Node_T::vec_dbl widths(3, 1.0);
  for (int i = 0; i < dim; i++)
    widths[i] = 10.0;

  // Create node
  Node_T::SP_node
    node(new Node_T(dim, "cartnode",
         Node_T::vec2_size_t(2*dim, Node_T::vec_size_t(dim-1, so)),  // space
         Node_T::vec_size_t(2*dim, po),                              // polar
         Node_T::vec_size_t(2*dim, ao),                              // azimuth
         Node_T::vec_size_t(2*dim, 0),                               // energy
         widths, db, mat, mesh));
  Ensure(node);
  return node;
}


} // end namespace erme_geometry

#endif // NODE_FIXTURE_HH_ 

//----------------------------------------------------------------------------//
//              end of file node_fixture.hh
//----------------------------------------------------------------------------//
