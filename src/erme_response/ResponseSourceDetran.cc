//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ResponseSourceDetran.cc
 *  @brief ResponseSourceDetran
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "ResponseSourceDetran.hh"
#include "ResponseSourceDetranDiffusion.t.hh"
#include "ResponseSourceDetran.t.hh"
#include "orthog/OrthogonalBasis.hh"

namespace erme_response
{

#define COUT(c) std::cout << c << std::endl;

using std::string;
using detran_orthog::OrthogonalBasis;

//----------------------------------------------------------------------------//
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

  std::cout << "********* BUILDING DETRAN<" << (int)D::dimension
            << "> SOURCE FOR NODE " << d_node->name() << std::endl;

  d_db = node->db();
  d_material = node->material();
  d_mesh = node->mesh();
  //d_db->display();
  // Ensure we compute boundary fluxes
  d_db->put<int>("compute_boundary_flux", 1);

  // We use a fixed boundary.
  d_db->put<std::string>("bc_west",     "fixed");
  d_db->put<std::string>("bc_east",     "fixed");
  d_db->put<std::string>("bc_south",    "fixed");
  d_db->put<std::string>("bc_north",    "fixed");
  d_db->put<std::string>("bc_bottom",   "fixed");
  d_db->put<std::string>("bc_top",      "fixed");

  // Create the solver and extract the boundary and quadrature
  d_solver = new Solver_T(d_db, d_material, d_mesh, true);
  d_solver->setup();       // Constructs quadrature, etc.
  d_solver->set_solver();  // Constructs the actual mg solver
  d_B = d_solver->boundary();
  d_quadrature = d_solver->quadrature();
//  d_quadrature->display();
//  THROW("ppp");
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

//----------------------------------------------------------------------------//
template <class B>
void ResponseSourceDetran<B>::
compute(SP_response response, const ResponseIndex &index)
{
//  std::cout << "COMPUTING RESPONSE FOR INDEX: " << index << std::endl;
  d_solver->boundary()->clear();
  d_solver->boundary()->clear_bc();
  set_boundary(index);
//  std::cout << "idx =" << index << std::endl;

//  std::cout << "********* OUTGOING BOUNDARY *********** " << std::endl;
//  d_B->display(false);
//  std::cout << "********* INCIDENT BOUNDARY *********** " << std::endl;
//  d_B->display(true);

  d_solver->solve(d_keff);
 // THROW("lala");
 // d_solver->state()->display();
//std::cout << "********* OUTGOING BOUNDARY *********** " << std::endl;
//  d_B->display(false);
//  std::cout << "********* INCIDENT BOUNDARY *********** " << std::endl;
//  d_B->display(true);

  expand(response, index);

  //response->display();
// THROW("lala");
}

//----------------------------------------------------------------------------//
template <class B>
void ResponseSourceDetran<B>::construct_basis()
{
  using namespace detran_orthog;
  using std::string;

  construct_energy_basis();

  //--------------------------------------------------------------------------//
  // SPACE
  //--------------------------------------------------------------------------//

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
        OrthogonalBasis::Parameters basis_s_p;
        basis_s_p.size  = d_mesh->number_cells(fd);
        basis_s_p.order = d_node->spatial_order(s, dim01);
        basis_s_p.orthonormal = true;
        d_basis_s[s][dim01] = OrthogonalBasis::Create(basis_s_type, basis_s_p);
      }
    }
  }

  // Skip angular stuff if diffusion.
  if (d_solver->discretization() == d_solver->DIFF) return;

  //--------------------------------------------------------------------------//
  // ANGLE
  //--------------------------------------------------------------------------//

  // Determine whether to expand angular flux or angular current
  if (d_db->check("erme_angular_expansion"))
    d_angular_flux = d_db->get<int>("erme_angular_expansion");

  if (B::D_T::dimension == 1)
  {
    construct_angular_basis_1D();
  }
  else if (B::D_T::dimension == 2)
  {
    construct_angular_basis_2D();
  }
  if (B::D_T::dimension == 3)
  {
    construct_angular_basis_3D();
  }
}

//----------------------------------------------------------------------------//
template <class B>
void ResponseSourceDetran<B>::construct_energy_basis()
{
  /*
   *  The default basis in energy is the discrete Diract delta basis, which
   *  is limited to full order expansions.  This is equivalent to a standard
   *  multigroup treatment.
   *
   *  As an alternative, any discrete basis can be used.  This includes
   *  the transformed DCT basis.
   *
   */
  //OrthogonalBasis::Factory_T::ShowRegistered();

  // default is a dirac basis, equivalent to full multigroup
  string basis_e_type = "ddf";
  if (d_db->check("basis_e_type"))
    basis_e_type = d_db->get<string>("basis_e_type");
  size_t ng = d_material->number_groups();
  size_t basis_e_order = ng - 1;
  if (d_db->check("basis_e_order"))
    basis_e_order = d_db->get<int>("basis_e_order");
  Insist(basis_e_order < ng, "Energy order must be <= number groups");
  OrthogonalBasis::Parameters basis_e_p;
  basis_e_p.order = basis_e_order;
  basis_e_p.size  = ng;
  basis_e_p.orthonormal = true;
  // transformed basis parameters
  if (basis_e_type == "trans")
  {
    Insist(d_db->check("basis_e_zeroth"),
      "Transformed energy basis requires the zeroth order be given.");
    basis_e_p.x = d_db->get<vec_dbl>("basis_e_zeroth");
    Insist(basis_e_p.x.size() == basis_e_p.size,
      "Inconsistent zeroth order size.");
    basis_e_p.transformed_key = "dlp";
    if (d_db->check("basis_e_transformed_option"))
    {
      basis_e_p.transformed_option =
        d_db->get<int>("basis_e_transformed_option");
    }
  }

  for (size_t s = 0; s < d_node->number_surfaces(); ++s)
    d_basis_e[s] = OrthogonalBasis::Create(basis_e_type, basis_e_p);

}

//----------------------------------------------------------------------------//
template <class B>
void ResponseSourceDetran<B>::construct_angular_basis_1D()
{
  string basis_p_type = "dlp";
  if (d_db->check("basis_p_type"))
    basis_p_type = d_db->get<string>("basis_p_type");
  COUT("BASIS P TYPE = " << basis_p_type)
  OrthogonalBasis::Parameters basis_p_p;
  basis_p_p.size = d_quadrature->number_angles_octant();
  basis_p_p.x    = d_quadrature->cosines(d_quadrature->MU);
  basis_p_p.qw   = d_quadrature->weights();
  basis_p_p.orthonormal = true;
  basis_p_p.lower_bound = 0.0;
  basis_p_p.upper_bound = 1.0;
  for (size_t s = 0; s < d_node->number_surfaces(); ++s)
  {
    basis_p_p.order = d_node->polar_order(s);
    d_basis_p[s] = OrthogonalBasis::Create(basis_p_type, basis_p_p);
  }
  COUT("lb=" << basis_p_p.lower_bound << " ub=" << basis_p_p.upper_bound);
  d_basis_p[0]->basis()->print_matlab("BP.out");
}

//----------------------------------------------------------------------------//
template <class B>
void ResponseSourceDetran<B>::construct_angular_basis_2D()
{
  /*
   *  Any possible combination of bases available in Detran is allowed.
   *  Using cheby for polar and/or clp for azimuth gets special treatment, since
   *  with proper transformation of variables, their combo yields aconservative
   *  basis as first presented by Rahnema and his students and examined more
   *  rigorously in my thesis.
   */

  SP_productquadrature q = d_quadrature;

  string basis_p_type = "dlp";
  if (d_db->check("basis_p_type"))
    basis_p_type = d_db->get<string>("basis_p_type");
  OrthogonalBasis::Parameters basis_p_p;
  basis_p_p.size = q->number_polar_octant();
  basis_p_p.orthonormal = true;
  basis_p_p.x.resize(basis_p_p.size, 0.0);
  basis_p_p.qw.resize(basis_p_p.size, 0.0);
  basis_p_p.lower_bound = 0.0;
  basis_p_p.upper_bound = 1.0;

  if (basis_p_type == "cheby")
  {
    // if cheby is used, we use the full -1 to 1 but skip odd orders (rather
    // define the basis from 0 to 1)
    basis_p_p.even_only   = true;
    basis_p_p.lower_bound = -1.0;
    basis_p_p.upper_bound =  1.0;
  }
  for (size_t p = 0; p < basis_p_p.size; ++p)
  {
    basis_p_p.x[p]  = q->cos_theta(p);
    basis_p_p.qw[p] = q->polar_weight(p);
  }

  string basis_a_type = "dlp";
  if (d_db->check("basis_p_type"))
    basis_a_type = d_db->get<string>("basis_a_type");
  OrthogonalBasis::Parameters basis_a_p;
  basis_a_p.size = 2 * q->number_azimuths_octant();
  basis_a_p.orthonormal = true;
  basis_a_p.x.resize(basis_a_p.size, 0.0);
  basis_a_p.qw.resize(basis_a_p.size, 0.0);
  if (basis_a_type == "clp")
  {
    // if clp is used for a, we transform via
    //     omega = cos(phi)
    // --> domega = -sin(phi) dphi
    // that negative is folded into the points, while the sin is factored into
    // the weights.
    for (size_t a = 0; a < basis_a_p.size / 2; ++a)
    {
      size_t b = basis_a_p.size - a - 1;
      basis_a_p.x[a]  = q->cos_phi(a);
      basis_a_p.x[b]  = -basis_a_p.x[a];
      basis_a_p.qw[a] = q->sin_phi(a) * q->azimuth_weight(a);
      basis_a_p.qw[b] = basis_a_p.qw[a];
    }
    basis_a_p.lower_bound = -1.0;
    basis_a_p.upper_bound =  1.0;
  }
  else
  {
    // otherwise, it's just a standard expansion over 0 to 2pi
    const double four_over_pi = 4.0 / detran_utilities::pi;
    for (size_t a = 0; a < basis_a_p.size / 2; ++a)
    {
      size_t b = basis_a_p.size - a - 1;
      basis_a_p.x[a]  = q->phi(a);
      basis_a_p.x[b]  = detran_utilities::pi - basis_a_p.x[a];
      basis_a_p.qw[a] = q->azimuth_weight(a) * four_over_pi;
      basis_a_p.qw[b] = basis_a_p.qw[a];
    }
    basis_a_p.lower_bound = 0.0;
    basis_a_p.upper_bound = detran_utilities::pi;
  }

  for (size_t s = 0; s < d_node->number_surfaces(); ++s)
  {
    basis_p_p.order = d_node->polar_order(s);
    basis_a_p.order = d_node->azimuthal_order(s);
    d_basis_p[s] = OrthogonalBasis::Create(basis_p_type, basis_p_p);
    d_basis_a[s] = OrthogonalBasis::Create(basis_a_type, basis_a_p);
  }
}

//----------------------------------------------------------------------------//
template <class B>
void ResponseSourceDetran<B>::construct_angular_basis_3D()
{
  /*
   *  For 3-D, we allow two options: DLP-DLP and a combination of
   *  cheby-clp + jacobi-clp.
   */

  SP_productquadrature q = d_quadrature;

  string basis_p_type = "dlp";
  if (d_db->check("basis_p_type"))
    basis_p_type = d_db->get<string>("basis_p_type");
  string basis_a_type = "dlp";
  if (d_db->check("basis_a_type"))
    basis_a_type = d_db->get<string>("basis_a_type");



  OrthogonalBasis::Parameters basis_p_p;
  basis_p_p.size = q->number_polar_octant();
  basis_p_p.orthonormal = true;
  basis_p_p.x.resize(basis_p_p.size, 0.0);
  basis_p_p.qw.resize(basis_p_p.size, 0.0);
  basis_p_p.lower_bound = 0.0;
  basis_p_p.upper_bound = 1.0;

}

//----------------------------------------------------------------------------//
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

//----------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//----------------------------------------------------------------------------//

template class ResponseSourceDetran<detran::BoundaryDiffusion<detran::_1D> >;
template class ResponseSourceDetran<detran::BoundaryDiffusion<detran::_2D> >;
template class ResponseSourceDetran<detran::BoundaryDiffusion<detran::_3D> >;
template class ResponseSourceDetran<detran::BoundarySN<detran::_1D> >;
template class ResponseSourceDetran<detran::BoundarySN<detran::_2D> >;
//template class ResponseSourceDetran<detran::BoundarySN<detran::_3D> >;
//template class ResponseSourceDetran<detran::BoundaryMOC<detran::_2D> >;

} // end namespace erme_response

//----------------------------------------------------------------------------//
//              end of file ResponseSourceDetran.cc
//----------------------------------------------------------------------------//
