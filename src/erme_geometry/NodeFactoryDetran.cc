//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   NodeFactoryDetran.cc
 * \brief  NodeFactoryDetran 
 * \author Jeremy Roberts
 * \date   Aug 23, 2012
 */
//---------------------------------------------------------------------------//

#include "NodeFactoryDetran.hh"
#include "CartesianNodeDetran.hh"

namespace erme_geometry
{

Node::SP_node NodeFactoryDetran::create_node(SP_db db,
                                             SP_material material,
                                             SP_mesh mesh)
{
  Require(db);
  Require(material);
  Require(mesh);

  typedef CartesianNodeDetran Node_T;

  // Get the dimension
  int dim = mesh->dimension();

  // Number of surfaces
  int ns  = 2 * dim;

  // Get the orders.  For now, assuming same order in all spatial directions.
  Insist(db->check("erme_spatial_order"), "Missing spatial order");
  Insist(db->check("erme_polar_order"), "Missing polar order");
  Insist(db->check("erme_azimuthal_order"), "Missing azimuthal order");
  // Nominally, we'll use delta's in energy (i.e. regular old multigroup).
  // However, one could also use expansions in a discrete basis.  This
  // would, in effect, allow mixed energy simulations.
  Insist(db->check("erme_energy_order"), "Missing energy order");
  int so = db->get<int>("erme_spatial_order");
  int po = db->get<int>("erme_polar_order");
  int ao = db->get<int>("erme_azimuthal_order");
  int eo = db->get<int>("erme_energy_order");
  Assert(so >= 0);
  Assert(po >= 0);
  Assert(ao >= 0);
  Assert(eo >= 0);
  Insist(eo < material->number_groups(),
         "The energy order must be smaller than the number of energy groups.");

  // Get the node name
  Insist(db->check("erme_node_name"), "Missing node name");
  std::string name = db->get<std::string>("erme_node_name");

  // Set the order vectors
  Node_T::vec2_size_t sov(ns, Node_T::vec_size_t(dim-1, so));
  Node_T::vec_size_t  pov(ns, po);
  Node_T::vec_size_t  aov(ns, ao);
  Node_T::vec_size_t  eov(ns, eo);

  // Origin
  Node_T::Point origin(0, 0, 0);

  // Width
  double w[] = {mesh->total_width_x(),
                mesh->total_width_y(),
                mesh->total_width_z()};
  Node_T::vec_dbl widths(3, 1.0);
  for (int i = 0; i < dim; i++)
    widths[i] = w[i];

  // Create node
  Node_T::SP_node node(new Node_T(dim, ns, 0, name, origin,
                                  sov, pov, aov, eov, widths,
                                  db, material, mesh));
  Ensure(node);
  return node;

}

} // end namespace erme_geometry

//---------------------------------------------------------------------------//
//              end of file NodeFactoryDetran.cc
//---------------------------------------------------------------------------//
