//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ResponseSourceDetranDiffusion.t.hh
 *  @brief ResponseSourceDetran diffusion specializations
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef erme_response_RESPONSESOURCEDETRANDIFFUSION_T_HH_
#define erme_response_RESPONSESOURCEDETRANDIFFUSION_T_HH_

#include "ResponseSourceDetran.hh"
#include "boundary/BoundaryDiffusion.hh"

namespace erme_response
{

//----------------------------------------------------------------------------//
// DIFFUSION SPECIALIZATIONS
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
template <>
void ResponseSourceDetran<detran::BoundaryDiffusion<detran::_1D> >::
set_boundary(const ResponseIndex &index)
{
  Boundary_T &B = *d_B;
  for (size_t g = 0; g < d_material->number_groups(); ++g)
  {
    Assert(d_basis_e[index.surface]);
    B(index.surface, g, B.IN) = (*d_basis_e[index.surface])(index.energy, g);
  }
  B.display(B.IN);
}

//----------------------------------------------------------------------------//
template <>
void ResponseSourceDetran<detran::BoundaryDiffusion<detran::_2D> >::
set_boundary(const ResponseIndex &index)
{
  Boundary_T &B = *d_B;
  double sign = 1.0;
  if ((index.surface == d_mesh->WEST   or
       index.surface == d_mesh->NORTH) and
      ((index.space0 + index.space1) % 2))
  {
    sign = -1.0;
  }
  size_t dim  = index.surface / 2;
  size_t dim0 = d_spatial_dim[dim][0];
  for (size_t g = 0; g < d_material->number_groups(); ++g)
  {
    double P_e = (*d_basis_e[index.surface])(index.energy, g);
    BoundaryTraits_T::value_type &b = B(index.surface, g, B.IN);
    for (size_t i = 0; i < d_mesh->number_cells(dim0); ++i)
    {
      double P_s0 = (*d_basis_s[index.surface][0])(index.space0, i);
      BoundaryValue_T::value(b, i, 0) =  sign * P_e * P_s0;
    }
  }
}

//----------------------------------------------------------------------------//
template <>
void ResponseSourceDetran<detran::BoundaryDiffusion<detran::_3D> >::
set_boundary(const ResponseIndex &index)
{
  Boundary_T &B = *d_B;
  double sign = 1.0;
  if ((index.surface == d_mesh->WEST  or
       index.surface == d_mesh->NORTH or
       index.surface == d_mesh->TOP)  and
      ((index.space0 + index.space1) % 2))
  {
    sign = -1.0;
  }
  for (size_t g = 0; g < d_material->number_groups(); ++g)
  {
    double P_e = (*d_basis_e[index.surface])(index.energy, g);
    size_t dim  = index.surface / 2;
    size_t dim0 = d_spatial_dim[dim][0];
    size_t dim1 = d_spatial_dim[dim][1];
    BoundaryTraits_T::value_type &b = B(index.surface, g, B.IN);
    for (size_t j = 0; j < d_mesh->number_cells(dim1); ++j)
    {
      double P_s1 = (*d_basis_s[index.surface][1])(index.space1, j);
      for (size_t i = 0; i < d_mesh->number_cells(dim0); ++i)
      {
        double P_s0 = (*d_basis_s[index.surface][0])(index.space0, i);
        BoundaryValue_T::value(b, i, j) = sign * P_e * P_s0 * P_s1;
      }
    }
  }
}

//----------------------------------------------------------------------------//
template <>
void ResponseSourceDetran<detran::BoundaryDiffusion<detran::_1D> >::
expand_boundary(SP_response          response,
                const ResponseIndex &index_i)
{
  const Boundary_T &B = *d_B;
  using namespace detran;
  typedef BoundaryTraits<_1D> B_T;
  typedef BoundaryValue<_1D> B_V;

  for (size_t surface = 0; surface < 2; ++surface)
  {
    //--------------------------------------------------------------------------//
    // CURRENT RESPONSE
    //--------------------------------------------------------------------------//

    vec_dbl b_e(d_material->number_groups(), 0.0);
    for (size_t g = 0; g < d_material->number_groups(); ++g)
    {
      b_e[g] = B(surface, g, B.OUT);
    }
    vec_dbl b_e2(d_basis_e[surface]->order() + 1, 0.0);
    d_basis_e[surface]->transform(b_e, b_e2);
    size_t nm = d_indexer->number_surface_moments(index_i.node, surface);
    for (size_t m = 0; m < nm; ++m)
    {
      ResponseIndex index_o = d_indexer->response_index(index_i.node, surface, m);
      response->boundary_response(index_o.nodal, index_i.nodal) = b_e[index_o.energy];
    }

    //------------------------------------------------------------------------//
    // LEAKAGE RESPONSE
    //------------------------------------------------------------------------//

    response->leakage_response(surface, index_i.nodal) = 0.0;
    for (size_t g = 0; g < d_material->number_groups(); ++g)
    {
      const B_T::value_type &bo = B(surface, g, B.OUT);
      const B_T::value_type &bi = B(surface, g, B.IN);
      response->leakage_response(surface, index_i.nodal) +=
        (B_V::value(bo, 0, 0) - B_V::value(bi, 0, 0));
    }
  }

}

//----------------------------------------------------------------------------//
template <>
void ResponseSourceDetran<detran::BoundaryDiffusion<detran::_2D> >::
expand_boundary(SP_response          response,
                const ResponseIndex &index_i)
{
  const Boundary_T &B = *d_B;
  using namespace detran;
  typedef BoundaryTraits<_2D> B_T;
  typedef BoundaryValue<_2D> B_V;

  for (size_t surface = 0; surface < 4; ++surface)
  {

    //------------------------------------------------------------------------//
    // CURRENT RESPONSE
    //------------------------------------------------------------------------//

    size_t o_s = d_basis_s[surface][0]->order();
    size_t o_e = d_basis_e[surface]->order();
    size_t n_g = d_material->number_groups();

    // Temporary response container
    vec2_dbl R(o_s + 1, vec_dbl(n_g, 0.0));
    vec_dbl Rs(o_s + 1, 0);

    // First expand in space, [spatial moments][energy groups]
    for (size_t g = 0; g < n_g; ++g)
    {
      const B_T::value_type &b = B(surface, g, B.OUT);
      d_basis_s[surface][0]->transform(b, Rs);
      for (size_t s = 0; s < Rs.size(); ++s)
      {
        R[s][g] = Rs[s];
      }
    }
    // Expand the result in energy, [spatial moments][energy moments]
    for (size_t s = 0; s <= o_s; ++s)
    {
      vec_dbl R_g(o_e + 1, 0.0);
      d_basis_e[surface]->transform(R[s], R_g);
      R[s] = R_g;
    }
    // Sign switch saves us from integrating in reverse direction.  This
    // assumes, of course, that basis functions are strictly even/odd.
    double sign = 1;
    if (surface == d_mesh->EAST or surface == d_mesh->SOUTH or
        surface == d_mesh->BOTTOM)
    {
      sign = -1;
    }
    // Fill the response container with only the *needed* values
    size_t nm = d_indexer->number_surface_moments(index_i.node, surface);
    for (size_t m = 0; m < nm; ++m)
    {
      ResponseIndex index_o = d_indexer->response_index(index_i.node, surface, m);
      double poo = std::pow(sign, index_o.space0);
      response->boundary_response(index_o.nodal, index_i.nodal) =
        poo * R[index_o.space0][index_o.energy];
    }

    //------------------------------------------------------------------------//
    // LEAKAGE RESPONSE
    //------------------------------------------------------------------------//

    size_t dim  = surface / 2;
    size_t dim0 = d_spatial_dim[dim][0];
    response->leakage_response(surface, index_i.nodal) = 0.0;
    for (size_t g = 0; g < d_material->number_groups(); ++g)
    {
      const B_T::value_type &bo = B(surface, g, B.OUT);
      const B_T::value_type &bi = B(surface, g, B.IN);

      for (size_t i = 0; i < bo.size(); ++i)
      {
        double dx = d_mesh->width(dim0, i);
        response->leakage_response(surface, index_i.nodal) +=
          dx * (B_V::value(bo, i) - B_V::value(bi, i));
      }
    }
  }

}

//----------------------------------------------------------------------------//
template <>
void ResponseSourceDetran<detran::BoundaryDiffusion<detran::_3D> >::
expand_boundary(SP_response          response,
                const ResponseIndex &index_i)
{
  const Boundary_T &B = *d_B;
  using namespace detran;
  typedef BoundaryTraits<_3D> B_T;
  typedef BoundaryValue<_3D> B_V;

  //--------------------------------------------------------------------------//
  // CURRENT RESPONSE
  //--------------------------------------------------------------------------//

  for (size_t surface = 0; surface < 6; ++surface)
  {
    size_t o_s0 = d_basis_s[surface][0]->order();
    size_t o_s1 = d_basis_s[surface][1]->order();
    size_t o_e  = d_basis_e[surface]->order();
    size_t n_g  = d_material->number_groups();

    size_t dim  = surface / 2;
    size_t dim0 = d_spatial_dim[dim][0];
    size_t dim1 = d_spatial_dim[dim][1];

    // Temporary response container
    vec3_dbl R(o_s0 + 1, vec2_dbl(o_s1 + 1, vec_dbl(n_g, 0.0)));

    for (size_t g = 0; g < n_g; ++g)
    {
      const B_T::value_type &b = B(surface, g, B.OUT);

      // Expand the first spatial variable
      vec2_dbl R_s0m_s1(o_s0 + 1, vec_dbl(b.size(), 0.0));
      vec_dbl s0m(o_s0 + 1, 0);
      for (size_t i = 0; i < b.size(); ++i) // For all s1
      {
        d_basis_s[surface][0]->transform(b[i], s0m);
        for (size_t j = 0; j < s0m.size(); ++j)
          R_s0m_s1[j][i] = s0m[j];
      }

      // Expand the second spatial variable
      vec_dbl s1m(o_s1 + 1, 0);
      for (size_t i = 0; i < R_s0m_s1.size(); ++i) // for all s0m
      {
        d_basis_s[surface][1]->transform(R_s0m_s1[i], s1m);
        for (size_t j = 0; j < s1m.size(); ++j)
          R[i][j][g] = s1m[j];
      }
    }

    // Expand the result in energy
    for (size_t i = 0; i < R.size(); ++i)
    {
      for (size_t j = 0; j < R[0].size(); ++j)
      {
        vec_dbl R_g(o_e + 1, 0.0);
        d_basis_e[surface]->transform(R[i][j], R_g);
        R[i][j] = R_g;
      }
    }

    // Sign switch saves us from integrating in reverse direction.  This
    // assumes, of course, that basis functions are strictly even/odd.
    double sign = spatial_sign(surface);

    // Fill the response container with only the *needed* values
    size_t nm = d_indexer->number_surface_moments(index_i.node, surface);
    for (size_t m = 0; m < nm; ++m)
    {
      ResponseIndex index_o =
        d_indexer->response_index(index_i.node, surface,m);
      double poo = std::pow(sign, index_o.space0 + index_o.space1);
      response->boundary_response(index_o.nodal, index_i.nodal) =
        poo * R[index_o.space0][index_o.space1][index_o.energy];
    }
  }

  //--------------------------------------------------------------------------//
  // LEAKAGE RESPONSE
  //--------------------------------------------------------------------------//

  for (size_t surface = 0; surface < 6; ++surface)
  {
    size_t dim  = surface / 2;
    size_t dim0 = d_spatial_dim[dim][0];
    size_t dim1 = d_spatial_dim[dim][1];
    response->leakage_response(surface, index_i.nodal) = 0.0;
    for (size_t g = 0; g < d_material->number_groups(); ++g)
    {
      const B_T::value_type &bo = B(surface, g, B.OUT);
      const B_T::value_type &bi = B(surface, g, B.IN);

      for (size_t j = 0; j < bo.size(); ++j)
      {
        double dx = d_mesh->width(dim1, j);
        for (size_t i = 0; i < bo[0].size(); ++i)
        {
          double dy = d_mesh->width(dim0, i);
          double da = dx * dy;
          response->leakage_response(surface, index_i.nodal) +=
            da * (B_V::value(bo, i, j) - B_V::value(bi, i, j));
        }
      }
    }
  }

}

} // end namespace erme_response

#endif /* erme_response_RESPONSESOURCEDETRANDIFFUSION_T_HH_ */

//----------------------------------------------------------------------------//
//              end of file ResponseSourceDetranDiffusion.t.hh
//----------------------------------------------------------------------------//
