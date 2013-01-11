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
template <class D>
ResponseSourceDetran<D>::ResponseSourceDetran(SP_node node)
  : ResponseSource(node)
  , d_basis_e(node->number_surfaces())
  , d_basis_s(node->number_surfaces(), vec_basis(D::dimension-1))
  , d_basis_a(node->number_surfaces())
  , d_basis_p(node->number_surfaces())
{
  // Preconditions
  Require(node->db());
  Require(node->material());
  Require(node->mesh());

  d_db = node->db();
  d_material = node->material();
  d_mesh = node->mesh();

  // Create the solver
  d_solver = new Solver_T(d_db, d_material, d_mesh, true);
  d_solver->setup();

  // Quadrature
  d_quadrature = d_solver->quadrature();

}

//---------------------------------------------------------------------------//
template <class D>
void ResponseSourceDetran<D>::construct_basis()
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

  size_t s = 0;
  for (size_t dim = 0; dim < D::dimension; ++dim)
  {
    for (size_t dir = 0; dir < 2; ++dir, ++s)
    {
      d_basis_s[s][dir] = new detran_orthog::
        DLP(d_node->spatial_order(s, dir), d_mesh->number_cells(dim));
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

  if (D::dimension > 1)
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
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

template class ResponseSourceDetran<detran::_1D>;
template class ResponseSourceDetran<detran::_2D>;
template class ResponseSourceDetran<detran::_3D>;

} // end namespace erme_response

//---------------------------------------------------------------------------//
//              end of file ResponseSourceDetran.cc
//---------------------------------------------------------------------------//
