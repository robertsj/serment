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
#include "postprocess/ReactionRates.hh"

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
  , d_expand_angular_flux(true)
  , d_spatial_dim(D::dimension, vec_size_t(D::dimension-1))
  , d_compute_nodal_power(false)
  , d_compute_pin_power(false)
{
  Require(node->db());
  Require(node->material());
  Require(node->mesh());

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

  std::string disc = "SN";
  if (d_solver->discretization() == d_solver->DIFF)
    disc = "DIFFUSION";
  else if (d_solver->discretization() == d_solver->MOC)
    disc = "MOC";

  std::cout << "********* BUILDING DETRAN<" << (int)D::dimension
            << "> SOURCE USING " << disc << " FOR NODE "
            << d_node->name() << std::endl;


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
   //if (index.surface < 2) return;

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
 // response->display();
// THROW("lala");
}

//----------------------------------------------------------------------------//
template <class B>
void ResponseSourceDetran<B>::construct_basis()
{
  using namespace detran_orthog;
  using std::string;

  //--------------------------------------------------------------------------//
  // ENERGY
  //--------------------------------------------------------------------------//

  construct_energy_basis();

  //--------------------------------------------------------------------------//
  // SPACE
  //--------------------------------------------------------------------------//

  string basis_s_type = "dlp";
  if (d_db->check("basis_s_type"))
    basis_s_type = d_db->get<string>("basis_s_type");
  //std::cout << " SPATIAL BASIS: " << basis_s_type << std::endl;

  // Loop over major dimension, direction (+/-), and secondary dimensions
  size_t s = 0;
  double w[3] = {d_mesh->total_width_x(),
                 d_mesh->total_width_y(),
                 d_mesh->total_width_z()};
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
        basis_s_p.x.resize(basis_s_p.size);
        basis_s_p.qw.resize(basis_s_p.size);
        double x = 0.0;
        double dx = 0.0;
        if (basis_s_type == "trans")
        {
          basis_s_p.transformed_key    = "dlp";
          basis_s_p.transformed_option = 1;
          for (int i = 0; i < basis_s_p.size; ++i)
          {
            basis_s_p.x[i]  = d_mesh->width(dim, i) / w[dim];
          }
        }
        else
        {
          for (int i = 0; i < basis_s_p.size; ++i)
          {
            basis_s_p.x[i]  = x + 0.5 * d_mesh->width(dim, i);
            basis_s_p.qw[i] = 2.0 * d_mesh->width(dim, i) / w[dim];
            x  += d_mesh->width(dim, i);
          }
        }
        basis_s_p.lower_bound = 0.0;
        basis_s_p.upper_bound = w[dim];
        d_basis_s[s][dim01] = OrthogonalBasis::Create(basis_s_type, basis_s_p);
      }
    }
  }

  // Skip angular stuff if diffusion.
  if (d_solver->discretization() == d_solver->DIFF)
  {
    return;
  }

  //--------------------------------------------------------------------------//
  // ANGLE
  //--------------------------------------------------------------------------//

  // Determine whether to expand angular flux or angular current
  if (d_db->check("erme_expand_angular_flux"))
  {
    d_expand_angular_flux = (1 == d_db->get<int>("erme_expand_angular_flux"));
  }

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
  d_basis_p[0]->basis()->print_matlab("BP.out");
}

//----------------------------------------------------------------------------//
template <class B>
void ResponseSourceDetran<B>::construct_angular_basis_2D()
{
  /*
   *  Any possible combination of bases available in Detran is
   *  allowed. Using cheby for polar and/or clp for azimuth gets
   *  special treatment, since with proper transformation of variables,
   *  their combo yields a conservative basis as first presented by
   *  Rahnema and his students and examined more rigorously in my thesis.
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
  basis_p_p.lower_bound = -1.0;
  basis_p_p.upper_bound =  1.0;

  if (basis_p_type == "cheby")
  {
    // if cheby is used, we use the full -1 to 1 but skip odd orders (rather
    // than define the basis from 0 to 1)
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
  if (d_db->check("basis_a_type"))
    basis_a_type = d_db->get<string>("basis_a_type");
  OrthogonalBasis::Parameters basis_a_p;
  basis_a_p.size = 2 * q->number_azimuths_octant();
  //basis_a_p.orthonormal = true;
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
    basis_a_p.orthonormal = true;
    basis_a_p.lower_bound = -1.0;
    basis_a_p.upper_bound =  1.0;
  }
  else
  {
    // otherwise, it's just a standard expansion over 0 to pi
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

  // For now, all surfaces are assigned the same basis
  for (size_t s = 0; s < d_node->number_surfaces(); ++s)
  {
    basis_p_p.order = d_node->polar_order(s);
    basis_a_p.order = d_node->azimuthal_order(s);
    d_basis_p[s] = OrthogonalBasis::Create(basis_p_type, basis_p_p);
    d_basis_a[s] = OrthogonalBasis::Create(basis_a_type, basis_a_p);
  }

  d_basis_p[0]->basis()->display(true);
  d_basis_p[0]->coefficients()->display();
  d_basis_p[0]->weights()->display();
}

//----------------------------------------------------------------------------//
template <class B>
void ResponseSourceDetran<B>::construct_angular_basis_3D()
{
  bool ORTHO = false;
  /*
   *  For 3-D, we allow two options: DLP-DLP or a combination of
   *  cheby-clp + jacobi-clp.
   */

  SP_productquadrature q = d_quadrature;

  string basis_a_type_v = "dlp";
  string basis_a_type_h = "dlp";
  string basis_p_type_v = "dlp";
  string basis_p_type_h = "dlp";

  if (d_db->check("basis_a_type"))
  {
    basis_a_type_v = d_db->get<string>("basis_a_type");
  }
  if (d_db->check("basis_a_type_h"))
  {
    basis_a_type_h = d_db->get<string>("basis_a_type_h");
  }
  if (d_db->check("basis_p_type"))
  {
    basis_p_type_v = d_db->get<string>("basis_p_type");
  }
  if (d_db->check("basis_p_type_h"))
  {
    basis_p_type_h = d_db->get<string>("basis_p_type_h");
  }

  Insist(basis_a_type_v == "dlp" || basis_a_type_v == "clp",
    "Only dlp or clp supported for 3D transport.");


  if (basis_a_type_v == "clp")
  {
    basis_p_type_v = "cheby";
  }
  if (basis_a_type_h == "clp")
  {
    //basis_p_type_h = "jacobi";
  }

  OrthogonalBasis::Parameters basis_p[2];
  basis_p[0].size =  2 * q->number_polar_octant();
  basis_p[1].size =  q->number_polar_octant();
  basis_p[0].lower_bound = -1.0;
  basis_p[1].lower_bound =  0.0;
  for (int i = 0; i < 2; ++i)
  {
    basis_p[i].orthonormal = ORTHO;
    basis_p[i].x.resize(basis_p[i].size, 0.0);
    basis_p[i].qw.resize(basis_p[i].size, 0.0);
    basis_p[i].upper_bound =  1.0;
  }

  for (size_t p = 0; p < q->number_polar_octant(); ++p)
  {
    size_t p0 = q->number_polar_octant() - p - 1;
    size_t p1 = q->number_polar_octant() + p    ;
    basis_p[0].x[p0]  = -q->cos_theta(p);
    basis_p[0].x[p1]  =  q->cos_theta(p);
    basis_p[0].qw[p0] =  2.0 * q->polar_weight(p) / detran_utilities::pi;
    basis_p[0].qw[p1] =  basis_p[0].qw[p0];
    //
    basis_p[1].x[p]   =  q->cos_theta(p);
    basis_p[1].qw[p]  =  q->polar_weight(p) * 2.0;
  }

  OrthogonalBasis::Parameters basis_a[2];
  basis_a[0].size = 2 * q->number_azimuths_octant();
  basis_a[1].size = 4 * q->number_azimuths_octant();
  for (int i = 0; i < 2; ++i)
  {
    basis_a[i].orthonormal = ORTHO;
    basis_a[i].x.resize(basis_a[i].size, 0.0);
    basis_a[i].qw.resize(basis_a[i].size, 0.0);
  }
  //basis_a[1].orthonormal = false;

  if (basis_a_type_v == "clp")
  {
    // if clp is used for a, we transform via
    //     omega = cos(phi)
    // --> domega = -sin(phi) dphi
    // that negative is folded into the points, while the sin is factored into
    // the weights.
    for (size_t a = 0; a < basis_a[0].size / 2; ++a)
    {
      size_t b = basis_a[0].size - a - 1;
      basis_a[0].x[a]  = q->cos_phi(a);
      basis_a[0].x[b]  = -basis_a[0].x[a];
      basis_a[0].qw[a] = 2*q->sin_phi(a) * q->azimuth_weight(a) / 2;
      basis_a[0].qw[b] = basis_a[0].qw[a];
    }
    basis_a[0].lower_bound = -1.0;
    basis_a[0].upper_bound =  1.0;
  }
  else
  {
    // otherwise, it's just a standard expansion over 0 to pi
    for (size_t a = 0; a < basis_a[0].size / 2; ++a)
    {
      size_t b = basis_a[0].size - a - 1;
      basis_a[0].x[a]  = q->phi(a);
      basis_a[0].x[b]  = detran_utilities::pi - basis_a[0].x[a];
      basis_a[0].qw[a] = 2 * q->azimuth_weight(a) / detran_utilities::pi;
      basis_a[0].qw[b] = basis_a[0].qw[a];
    }
    basis_a[0].lower_bound = 0.0;
    basis_a[0].upper_bound = detran_utilities::pi;
  }

  // For horizontal surfaces, CLP is done over 2pi
  for (size_t a = 0; a < basis_a[1].size / 4; ++a)
  {
    size_t b = basis_a[1].size / 2 - a - 1;
    size_t c = basis_a[1].size / 2 + a;
    size_t d = basis_a[1].size - a - 1;
    basis_a[1].x[a]  = q->phi(a);
    basis_a[1].x[b]  = detran_utilities::pi - basis_a[1].x[a];
    basis_a[1].x[c]  = detran_utilities::pi + basis_a[1].x[a];
    basis_a[1].x[d]  = 2.0 * detran_utilities::pi - basis_a[1].x[a];
    basis_a[1].qw[a] = 2 * q->azimuth_weight(a) / (2 * detran_utilities::pi);
    basis_a[1].qw[b] = basis_a[1].qw[a];
    basis_a[1].qw[c] = basis_a[1].qw[a];
    basis_a[1].qw[d] = basis_a[1].qw[a];
  }
  basis_a[1].lower_bound = 0.0;
  basis_a[1].upper_bound = 2.0 * detran_utilities::pi;

  // Vertical surfaces
  for (size_t s = 0; s < 4; ++s)
  {
    basis_p[0].order = d_node->polar_order(s);
    basis_a[0].order = d_node->azimuthal_order(s);
    d_basis_p[s] = OrthogonalBasis::Create(basis_p_type_v, basis_p[0]);
    d_basis_a[s] = OrthogonalBasis::Create(basis_a_type_v, basis_a[0]);
  }

  // Horizontal surfaces
  for (size_t s = 4; s < 6; ++s)
  {
    basis_p[1].order = d_node->polar_order(s);
    basis_a[1].order = d_node->azimuthal_order(s);
    d_basis_p[s] = OrthogonalBasis::Create(basis_p_type_h, basis_p[1]);
    d_basis_a[s] = OrthogonalBasis::Create(basis_a_type_h, basis_a[1]);
//    COUT("AZIMUTH")
//    for (int i = 0; i < basis_a[1].x.size(); ++i)
//    {
//      COUT(s << " " <<  i << " " << basis_a[1].x[i] << " " <<  basis_a[1].qw[i])
//    }
//    COUT("POLAR")
//    for (int i = 0; i < basis_p[1].x.size(); ++i)
//    {
//      COUT(s << " " << i << " " << basis_p[1].x[i] << " " << basis_p[1].qw[i])
//    }
  }
//
//  COUT("LALALALA")


    // Debug check on the quadratures to ensure they integrate the
    // partial current the way I expect.
  if (1)
  {
    using detran_utilities::soft_equiv;
    vec_dbl blah(basis_p[0].size, 1.0);
    vec_dbl blah_t(basis_p[0].order + 1, 1.0);
    d_basis_p[0]->transform(blah, blah_t);
    std::cout << " --> " << blah_t[0] << std::endl;
    //Assert(soft_equiv(blah_t[0], detran_utilities::pi / 2.0, 1e-9));

    blah.resize(basis_p[1].size, 1.0);
    blah_t.resize(basis_p[1].order + 1, 1.0);
    d_basis_p[4]->transform(blah, blah_t);
    std::cout << " --> " << blah_t[0] << std::endl;
    //Assert(soft_equiv(blah_t[0], 0.5, 1e-9));

    blah.resize(basis_a[0].size, 1.0);
    blah_t.resize(basis_a[0].order + 1, 1.0);
    d_basis_a[0]->transform(blah, blah_t);
    std::cout << " --> " << blah_t[0] << std::endl;

    blah.resize(basis_a[1].size, 1.0);
    blah_t.resize(basis_a[1].order + 1, 1.0);
    d_basis_a[4]->transform(blah, blah_t);
    std::cout << " --> " << blah_t[0] << std::endl;

  }
//  std::cout << " BASIS A " << std::endl;
//  d_basis_a[4]->basis()->display(true);
//  d_basis_a[5]->basis()->display(true);
//  THROW("lala");
}

//----------------------------------------------------------------------------//
template <class B>
void ResponseSourceDetran<B>::expand(SP_response          response,
                                     const ResponseIndex &index_i)
{
  typename Solver_T::SP_state state = d_solver->state();


  // Fission and absorption responses for the balance problem
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

  // Other volume responses
  detran_postprocess::ReactionRates
    rates(d_material, d_mesh, d_solver->state());
  vec_dbl tmp = rates.region_power("NODAL", -1.0);
  Assert(tmp.size() == 1);
  response->nodal_power(index_i.nodal) = tmp[0];
 // std::cout << " index=" << index_i.nodal << " nodal power = " << tmp[0] << std::endl;
  if (d_node->number_pins() > 0)
  {
    vec_dbl ppwr = rates.region_power("PINS", -1.0);
    Assert(ppwr.size() == response->number_pins());
    for (size_t p = 0; p < response->number_pins(); ++p)
      response->pin_power(p, index_i.nodal) = ppwr[p];
  }

  // Surface responses
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
template class ResponseSourceDetran<detran::BoundarySN<detran::_3D> >;
//template class ResponseSourceDetran<detran::BoundaryMOC<detran::_2D> >;

} // end namespace erme_response

//----------------------------------------------------------------------------//
//              end of file ResponseSourceDetran.cc
//----------------------------------------------------------------------------//
