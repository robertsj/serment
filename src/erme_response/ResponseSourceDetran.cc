//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ResponseSourceDetran.cc
 *  @brief  ResponseSourceDetran
 *  @author Jeremy Roberts
 *  @date   Jan 9, 2013
 */
//---------------------------------------------------------------------------//

#include "ResponseSourceDetran.hh"
#include "ResponseSourceDetranDiffusion.t.hh"
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
  , d_angular_flux(true)
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

  // We use a fixed boundary.
//  d_db->put<std::string>("bc_west",     "fixed");
//  d_db->put<std::string>("bc_east",     "fixed");
//  d_db->put<std::string>("bc_south",    "fixed");
//  d_db->put<std::string>("bc_north",    "fixed");
//  d_db->put<std::string>("bc_bottom",   "fixed");
//  d_db->put<std::string>("bc_top",      "fixed");

  // Create the solver and extract the boundary and quadrature
  d_solver = new Solver_T(d_db, d_material, d_mesh, true);
  d_solver->setup();       // Constructs quadrature, etc.
  d_solver->set_solver();  // Constructs the actual mg solver
  d_B = d_solver->boundary();
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
  //std::cout << "COMPUTING RESPONSE FOR INDEX: " << index << std::endl;
  d_solver->boundary()->clear();
  d_solver->boundary()->clear_bc();
  set_boundary(index);
  //THROW("lala");
//  std::cout << "********* OUTGOING BOUNDARY *********** " << std::endl;
//  d_B->display(false);
//  std::cout << "********* INCIDENT BOUNDARY *********** " << std::endl;
//  d_B->display(true);
  d_solver->solve(d_keff);
//  d_solver->state()->display();
//  std::cout << "********* OUTGOING BOUNDARY *********** " << std::endl;
//  d_B->display(false);
//  std::cout << "********* INCIDENT BOUNDARY *********** " << std::endl;
//  d_B->display(true);

  expand(response, index);
//  response->display();
//  THROW("lala");
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
      }
    }
  }

  // Skip angular stuff if diffusion.
  if (d_solver->discretization() == d_solver->DIFF) return;

  //-------------------------------------------------------------------------//
  // ANGLE
  //-------------------------------------------------------------------------//

  /*
   * Here, we need to do a bit more logic. If we do the DLP approach, then
   * a split of the polar and azimuth is useful.  In 1-D, it's only polar.
   * In 2-D, the split is the same on all sides (we just need to iterate
   * over the correct azimuths).  In 3-D, things change a bit.  On the
   * vertical surfaces, we have a half-azimuthal range and a full
   * polar range.  On the horizontal surfaces, it's a full azimuthal
   * range and a half polar range.  Note, in this way, the "polar" vs
   * "azimuthal" distinction is with respect to the *actual* geometry,
   * not with respect to the incident surface coordinates.
   *
   * If we go the route of continuous polynomials, it seems much more
   * straightforward to employ the local coordinate system.  Then, the
   * polar angle is with respect to the incident normal.  Then in 1-D,
   * it's just the Jacobi polynomial.  In 2-D/3-D, it's the Jacobi
   * polynomial and Legendre in the azimuth.  As long as we use a
   * quadrature that integrates these exactly, we should be fine.
   * It might be necessary to define a product quad with Legendre
   * in the azimuth, or we could just use level symmetric.  My guess is
   * the accuracy to integrate the Jacobi/Legendre moments is similar
   * to that needed for spherical harmonics.
   *
   */

  // Determine whether to expand angular flux or angular current
  if (d_db->check("erme_angular_expansion"))
    d_angular_flux = d_db->get<int>("erme_angular_expansion");

  // Polar
  string basis_p_type = "dlp";
  if (d_db->check("basis_p_type"))
    basis_p_type = d_db->get<string>("basis_p_type");
  //std::cout << " POLAR TYPE = " << basis_p_type << std::endl;

  if (B::D_T::dimension == 1)
  {
    size_t np = d_quadrature->number_angles_octant();
    for (size_t s = 0; s < d_node->number_surfaces(); ++s)
    {
      if (basis_p_type == "dlp")
      {
        d_basis_p[s] = new detran_orthog::
          DLP(d_node->polar_order(s), np, true);
      }
      else if (basis_p_type == "jacobi")
      {
        vec_dbl mu = d_quadrature->cosines(d_quadrature->MU);
        vec_dbl wt = d_quadrature->weights();
        d_basis_p[s] = new detran_orthog::
          Jacobi01(d_node->polar_order(s), mu, wt, 0.0, 1.0);
      }
      else if (basis_p_type == "legendre")
      {
        vec_dbl mu = d_quadrature->cosines(d_quadrature->MU);
        vec_dbl wt = d_quadrature->weights();
        d_basis_p[s] = new detran_orthog::
          CLP(d_node->polar_order(s), mu, wt, 0.0, 1.0);
      }
    }
  }
  else if (B::D_T::dimension == 2)
  {
    // Azimuth
    string basis_a_type = "dlp";
    if (d_db->check("basis_a_type"))
      basis_a_type = d_db->get<string>("basis_a_type");

    if (basis_p_type == "dlp")
    {
      // Use DLP for the *physical* polar and azimuth.
      SP_productquadrature q = d_quadrature;
      size_t np = q->number_polar_octant();
      size_t na = 2 * q->number_azimuths_octant();
      for (size_t s = 0; s < d_node->number_surfaces(); ++s)
      {
        d_basis_p[s] = new detran_orthog::DLP(d_node->polar_order(s),     np, true);
        d_basis_a[s] = new detran_orthog::DLP(d_node->azimuthal_order(s), na, true);
      }
    }
    else if (basis_p_type == "legendre")
    {
      // Use continuous Legendre for the polar and azimuth.
      SP_productquadrature q = d_quadrature;
      size_t np = q->number_polar_octant();
      size_t na = 2 * q->number_azimuths_octant();
      vec_dbl xi(np, 0);
      vec_dbl w_p(np, 0);
      for (size_t p = 0; p < np; ++p)
      {
        xi[p]  = q->cos_theta(p);
        w_p[p] = q->polar_weight(p);
      }
      vec_dbl phi(na, 0);
      vec_dbl w_a(na, 0);
      for (size_t a = 0; a < na/2; ++a)
      {
        phi[a]      = q->phi(a);
        phi[na-a-1] = detran_utilities::pi - q->phi(a);
        w_a[a]      = q->azimuth_weight(a) * 2.0 /  1.570796326794897;;
        w_a[na-a-1] = q->azimuth_weight(a) * 2.0 /  1.570796326794897;
      }
      for (size_t s = 0; s < d_node->number_surfaces(); ++s)
      {
        d_basis_p[s] = new detran_orthog::
          CLP(d_node->polar_order(s), xi, w_p, 0.0, 1.0);
        d_basis_a[s] = new detran_orthog::
          CLP(d_node->azimuthal_order(s), phi, w_a, 0.0, detran_utilities::pi);
      }
    }
    else if (basis_p_type == "cheby")
    {
      SP_productquadrature q = d_quadrature;
      size_t np = q->number_polar_octant();
      size_t na = 2 * q->number_azimuths_octant();
      vec_dbl xi(np, 0);
      vec_dbl w_p(np, 0);
      for (size_t p = 0; p < np; ++p)
      {
        xi[p]  = q->cos_theta(p);
        w_p[p] = q->polar_weight(p);
        //std::cout << " p = " << p << " w = " << w_p[p] << std::endl;
      }
      vec_dbl cos_phi(na, 0);
      vec_dbl sin_phi(na, 0);
      vec_dbl w_a(na, 0);
      for (size_t a = 0; a < na/2; ++a)
      {
        cos_phi[a]      = q->cos_phi(a);
        cos_phi[na-a-1] = -q->cos_phi(a);
        sin_phi[a]      = q->sin_phi(a);
        sin_phi[na-a-1] = q->sin_phi(a);
        //std::cout << " a = " << a << " cos_phi = " << q->cos_phi(a) << std::endl;
        w_a[a]      = q->sin_phi(a) * q->azimuth_weight(a);// /  1.570796326794897;
        w_a[na-a-1] = q->sin_phi(a) * q->azimuth_weight(a);// /  1.570796326794897;
      }
      for (size_t s = 0; s < d_node->number_surfaces(); ++s)
      {
        d_basis_p[s] = new detran_orthog::
          ChebyshevU(d_node->polar_order(s), xi, w_p, -1.0, 1.0, true);
        d_basis_a[s] = new detran_orthog::
          CLP(d_node->azimuthal_order(s), cos_phi, w_a, -1.0, 1.0);
      }
//      d_basis_p[0]->weights()->display();
//      d_basis_p[0]->basis()->display();
//      d_basis_a[0]->weights()->display();
//      d_basis_a[0]->basis()->display();

    }
    else if (basis_p_type == "jacobi")
    {
      // Use Jacobi for the polar w/r to the incident.
      // Note, this one

      size_t axis = s / 2;
      vec_dbl mu = d_quadrature->cosines(axis);
      vec_dbl wt = d_quadrature->weights();
      d_basis_p[s] = new detran_orthog::
        Jacobi01(d_node->polar_order(s), mu, wt, 0.0, 1.0);
    }
    else
    {
      THROW("INVALID BASIS");
    }

  }
  if (B::D_T::dimension == 3)
  {
    // Azimuth
    string basis_a_type = "dlp";
    if (d_db->check("basis_a_type"))
      basis_a_type = d_db->get<string>("basis_a_type");
    for (size_t s = 0; s < d_node->number_surfaces(); ++s)
      d_basis_a[s] = new detran_orthog::DLP(d_node->azimuthal_order(s), 0);
  }

}

//---------------------------------------------------------------------------//
template <class B>
void ResponseSourceDetran<B>::expand(SP_response response,
                                     const ResponseIndex &index_i)
{
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
  expand_boundary(response, index_i);
}

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

template class ResponseSourceDetran<detran::BoundaryDiffusion<detran::_1D> >;
template class ResponseSourceDetran<detran::BoundaryDiffusion<detran::_2D> >;
template class ResponseSourceDetran<detran::BoundaryDiffusion<detran::_3D> >;
template class ResponseSourceDetran<detran::BoundarySN<detran::_1D> >;
template class ResponseSourceDetran<detran::BoundarySN<detran::_2D> >;
//template class ResponseSourceDetran<detran::BoundarySN<detran::_3D> >;
//template class ResponseSourceDetran<detran::BoundaryMOC<detran::_2D> >;

} // end namespace erme_response

//---------------------------------------------------------------------------//
//              end of file ResponseSourceDetran.cc
//---------------------------------------------------------------------------//
