//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ResponseSourceDetran.cc
 *  @brief  ResponseSourceDetran
 *  @author Jeremy Roberts
 *  @date   Jan 9, 2013
 */
//---------------------------------------------------------------------------//

#include "ResponseSourceDetran.hh"
#include "ResponseSourceDetran.t.hh"
#include "orthog/detran_orthog.hh"

namespace erme_response
{

//---------------------------------------------------------------------------//
template <class B>
ResponseSourceDetran<B>::ResponseSourceDetran(SP_node node,
                                              SP_indexer indexer)
  : ResponseSource(node, indexer)
  , d_basis_e(node->number_surfaces())
  , d_basis_s(node->number_surfaces(), vec_basis(D::dimension-1))
  , d_basis_a(node->number_surfaces())
  , d_basis_p(node->number_surfaces())
  , d_spatial_dim(D::dimension, vec_size_t(D::dimension-1))
{
  Require(node->db());
  Require(node->material());
  Require(node->mesh());

  std::cout << "********* BUILDING DETRAN<" << D::dimension
            << "> SOURCE FOR NODE " << d_node->name() << std::endl;

  d_db = node->db();
  d_material = node->material();
  d_mesh = node->mesh();

  // Ensure we compute boundary fluxes
  d_db->put<int>("compute_boundary_flux", 1);

  // Create the solver
  d_solver = new Solver_T(d_db, d_material, d_mesh, true);
  d_solver->setup();       // Constructs quadrature, etc.
  d_solver->set_solver();  // Constructs the actual mg solver

  // Quadrature.  Remember, this may be NULL.
  d_quadrature = d_solver->quadrature();

  // Spatial dimensions in play.  For example, when expanding
  // on an x-directed surface, y and z are in play.
  if (D::dimension == 2)
  {
    d_spatial_dim[0][0] = 1;
    d_spatial_dim[1][0] = 0;
  }
  else if (D::dimension == 3)
  {
    d_spatial_dim[0][0] = 1;
    d_spatial_dim[0][1] = 2;
    d_spatial_dim[1][0] = 0;
    d_spatial_dim[1][1] = 2;
    d_spatial_dim[2][0] = 0;
    d_spatial_dim[2][1] = 1;
  }

  construct_basis();
}

//---------------------------------------------------------------------------//
template <class B>
void ResponseSourceDetran<B>::
compute(SP_response response, const ResponseIndex &index)
{
  using namespace detran;

  //std::cout << "COMPUTING RESPONSE FOR INDEX: " << index << std::endl;

  //-------------------------------------------------------------------------//
  // SET INCIDENT CONDITION
  //-------------------------------------------------------------------------//

  d_solver->boundary()->clear();
  typename Boundary_T::SP_boundary b = d_solver->boundary();
  set_boundary(*b, index);

  //-------------------------------------------------------------------------//
  // SOLVE THE RESPONSE EQUATION
  //-------------------------------------------------------------------------//

  d_solver->solve(d_keff);

  //-------------------------------------------------------------------------//
  // EXPAND THE RESPONSE
  //-------------------------------------------------------------------------//

  expand(*b, response, index);

}

//---------------------------------------------------------------------------//
template <class B>
void ResponseSourceDetran<B>::construct_basis()
{
  // @todo MOC is not incorporated yet

  using std::string;

  //-------------------------------------------------------------------------//
  // ENERGY
  //-------------------------------------------------------------------------//

  string basis_e_type = "ddf";
  if (d_db->check("basis_e_type"))
    basis_e_type = d_db->get<string>("basis_e_type");
  size_t ng = d_material->number_groups();
  for (size_t s = 0; s < d_node->number_surfaces(); ++s)
    d_basis_e[s] = new detran_orthog::DDF(d_node->energy_order(s), ng);

  //-------------------------------------------------------------------------//
  // SPACE
  //-------------------------------------------------------------------------//

  string basis_s_type = "dlp";
  if (d_db->check("basis_s_type"))
    basis_s_type = d_db->get<string>("basis_s_type");

  // Loop over major dimension, direction (+/-), and secondary dimensions
  size_t s = 0;
  for (size_t dim = 0; dim < D::dimension; ++dim) // INCIDENT DIMENSION
  {
    for (size_t dir = 0; dir < 2; ++dir, ++s)     // LEFT or RIGHT
    {
      for (size_t dim01 = 0; dim01 < D::dimension - 1; ++dim01) // FREE DIMENSIONS
      {
        size_t fd = d_spatial_dim[dim][dim01];
        d_basis_s[s][dim01] = new detran_orthog::
          DLP(d_node->spatial_order(s, dim01),
              d_mesh->number_cells(fd), true);
        std::cout << " s=" << s
                  << " dim = " << dim
                  << " dim01=" << dim01
                  << " fd = " << fd
                  << std::endl;
        d_basis_s[s][dim01]->basis()->display();
      }
    }
  }

  // Skip angular stuff if diffusion.
  if (d_solver->discretization() == d_solver->DIFF) return;

  //-------------------------------------------------------------------------//
  // ANGLE
  //-------------------------------------------------------------------------//

  size_t np = 0;
  size_t na = 0;

  if (D::dimension == 1)
  {
    np = d_quadrature->number_angles_octant();
  }
  else
  {
    SP_productquadrature q = d_quadrature;
    np = q->number_polar_octant();
    na = 2 * q->number_azimuths_octant();
  }

  // Polar
  string basis_p_type = "dlp";
  if (d_db->check("basis_p_type"))
    basis_p_type = d_db->get<string>("basis_p_type");
  for (size_t s = 0; s < d_node->number_surfaces(); ++s)
    d_basis_p[s] = new detran_orthog::DLP(d_node->polar_order(s), np);

  if (B::D_T::dimension > 1)
  {
    // Azimuth
    string basis_a_type = "dlp";
    if (d_db->check("basis_a_type"))
      basis_a_type = d_db->get<string>("basis_a_type");
    for (size_t s = 0; s < d_node->number_surfaces(); ++s)
      d_basis_a[s] = new detran_orthog::DLP(d_node->azimuthal_order(s), na);
  }

}

//---------------------------------------------------------------------------//
template <class B>
void ResponseSourceDetran<B>::expand_flux(SP_response           response,
                                          const ResponseIndex  &index_i)
{
  //-------------------------------------------------------------------------//
  // FISSION & ABSORPTION
  //-------------------------------------------------------------------------//

  typename Solver_T::SP_state state = d_solver->state();
  const vec_int &mat_map = d_mesh->mesh_map("MATERIAL");
  response->fission_response(index_i.nodal) = 0.0;
  response->absorption_response(index_i.nodal) = 0.0;
  for (size_t g = 0; g < d_material->number_groups(); ++g)
  {
    for (size_t i = 0; i < d_mesh->number_cells(); ++i)
    {
      double phi_times_volume = d_mesh->volume(i) * state->phi(g)[i];
      response->fission_response(index_i.nodal) +=
         phi_times_volume * d_material->nu_sigma_f(mat_map[i], g);
      response->absorption_response(index_i.nodal) +=
         phi_times_volume * d_material->sigma_a(mat_map[i], g);
    }
  }
}

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

template class ResponseSourceDetran<detran::BoundaryDiffusion<detran::_1D> >;
template class ResponseSourceDetran<detran::BoundaryDiffusion<detran::_2D> >;
template class ResponseSourceDetran<detran::BoundaryDiffusion<detran::_3D> >;
template class ResponseSourceDetran<detran::BoundarySN<detran::_1D> >;
//template class ResponseSourceDetran<detran::BoundarySN<detran::_2D> >;
//template class ResponseSourceDetran<detran::BoundarySN<detran::_3D> >;
//template class ResponseSourceDetran<detran::BoundaryMOC<detran::_2D> >;

} // end namespace erme_response

//---------------------------------------------------------------------------//
//              end of file ResponseSourceDetran.cc
//---------------------------------------------------------------------------//
