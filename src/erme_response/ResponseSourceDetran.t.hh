//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ResponseSourceDetran.t.hh
 *  @brief  ResponseSourceDetran
 *  @author Jeremy Roberts
 *  @date   Jan 9, 2013
 */
//---------------------------------------------------------------------------//

#ifndef erme_response_RESPONSESOURCEDETRAN_T_HH_
#define erme_response_RESPONSESOURCEDETRAN_T_HH_

#include "ResponseSourceDetran.hh"
#include "boundary/BoundaryDiffusion.hh"
#include "boundary/BoundarySN.hh"
#include "boundary/BoundaryMOC.hh"
#include "boundary/BoundaryTraits.hh"

namespace erme_response
{

// \todo Need to consolidate code where possible

//---------------------------------------------------------------------------//
template <class B>
void ResponseSourceDetran<B>::set_boundary(B& boundary,
                                           const ResponseIndex &index)
{
  THROW("NOT IMPLEMENTED");
}
//---------------------------------------------------------------------------//
template <class B>
void ResponseSourceDetran<B>::expand(const B             &boundary,
                                     SP_response          response,
                                     const ResponseIndex &index_i)
{
  THROW("NOT IMPLEMENTED");
}

//---------------------------------------------------------------------------//
// DIFFUSION SPECIALIZATIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
template <>
void ResponseSourceDetran<detran::BoundaryDiffusion<detran::_1D> >::
set_boundary(Boundary_T &boundary, const ResponseIndex &index)
{
  for (size_t g = 0; g < d_material->number_groups(); ++g)
  {
    BoundaryTraits_T::value_type &b = boundary(index.surface, g, boundary.IN);
    Assert(d_basis_e[index.surface]);
    BoundaryValue_T::value(b, 0, 0) = (*d_basis_e[index.surface])(index.energy, g);
  }
}

//---------------------------------------------------------------------------//
template <>
void ResponseSourceDetran<detran::BoundaryDiffusion<detran::_2D> >::
set_boundary(Boundary_T &boundary, const ResponseIndex &index)
{
  using namespace detran;
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
    BoundaryTraits<_2D>::value_type
      &b = boundary(index.surface, g, boundary.IN);
    for (size_t i = 0; i < d_mesh->number_cells(dim0); ++i)
    {
      double P_s0 = (*d_basis_s[index.surface][0])(index.space0, i);
      BoundaryValue<_2D>::value(b, i, 0) =  sign * P_e * P_s0;
    }
  }
}

//---------------------------------------------------------------------------//
template <>
void ResponseSourceDetran<detran::BoundaryDiffusion<detran::_3D> >::
set_boundary(Boundary_T &boundary, const ResponseIndex &index)
{
  using namespace detran;
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
    BoundaryTraits<_3D>::value_type
      &b = boundary(index.surface, g, boundary.IN);
//    std::cout << " 0=" << b.size()
//              << " 1=" << b[0].size()
//              << " dim =" << dim
//              << " dim0=" << dim0
//              << " dim1=" << dim1
//              << std::endl;
    for (size_t j = 0; j < d_mesh->number_cells(dim1); ++j)
    {
      double P_s1 = (*d_basis_s[index.surface][1])(index.space1, j);
      for (size_t i = 0; i < d_mesh->number_cells(dim0); ++i)
      {
        double P_s0 = (*d_basis_s[index.surface][0])(index.space0, i);

        BoundaryValue<_3D>::value(b, i, j) = sign * P_e * P_s0 * P_s1;
      }
    }
  }
  //std::cout << " VALUE SET " << std::endl;
}

//---------------------------------------------------------------------------//
template <>
void ResponseSourceDetran<detran::BoundaryDiffusion<detran::_1D> >::
expand(const Boundary_T    &boundary,
       SP_response          response,
       const ResponseIndex &index_i)
{
  using namespace detran;
  typedef BoundaryTraits<_1D> B_T;
  typedef BoundaryValue<_1D> B_V;

  //-------------------------------------------------------------------------//
  // CURRENT RESPONSE
  //-------------------------------------------------------------------------//

  for (size_t surface = 0; surface < 2; ++surface)
  {
    // Temporary energy-dependent vector
    vec_dbl b_e(d_material->number_groups(), 0.0);

    for (size_t g = 0; g < d_material->number_groups(); ++g)
    {
      const B_T::value_type &b = boundary(surface, g, boundary.OUT);
      b_e[g] = B_V::value(b, 0, 0);
    }
    vec_dbl b_e2(d_basis_e[surface]->order() + 1, 0.0);
    d_basis_e[surface]->transform(b_e, b_e2);
    size_t nm = d_indexer->number_surface_moments(index_i.node, surface);
    for (size_t m = 0; m < nm; ++m)
    {
      ResponseIndex index_o = d_indexer->response_index(index_i.node, surface, m);
      response->boundary_response(index_o.nodal, index_i.nodal) = b_e[index_o.energy];
    }
  }

  //-------------------------------------------------------------------------//
  // LEAKAGE RESPONSE
  //-------------------------------------------------------------------------//

  for (size_t surface = 0; surface < 2; ++surface)
  {
    response->leakage_response(surface, index_i.nodal) = 0.0;
    for (size_t g = 0; g < d_material->number_groups(); ++g)
    {
      const B_T::value_type &bo = boundary(surface, g, boundary.OUT);
      const B_T::value_type &bi = boundary(surface, g, boundary.IN);
      response->leakage_response(surface, index_i.nodal) +=
        (B_V::value(bo, 0, 0) - B_V::value(bi, 0, 0));
    }
  }

  //-------------------------------------------------------------------------//
  // FLUX RESPONSES
  //-------------------------------------------------------------------------//

  expand_flux(response, index_i);

}

//---------------------------------------------------------------------------//
template <>
void ResponseSourceDetran<detran::BoundaryDiffusion<detran::_2D> >::
expand(const Boundary_T    &boundary,
       SP_response          response,
       const ResponseIndex &index_i)
{
  using namespace detran;
  typedef BoundaryTraits<_2D> B_T;
  typedef BoundaryValue<_2D> B_V;

  //-------------------------------------------------------------------------//
  // CURRENT RESPONSE
  //-------------------------------------------------------------------------//

  for (size_t surface = 0; surface < 4; ++surface)
  {
    size_t o_s = d_basis_s[surface][0]->order();
    size_t o_e = d_basis_e[surface]->order();
    size_t n_g = d_material->number_groups();

    // Temporary response container
    vec2_dbl R(o_s + 1, vec_dbl(n_g, 0.0));
    vec_dbl Rs(o_s + 1, 0);

    // First expand in space, [spatial moments][energy groups]
    for (size_t g = 0; g < n_g; ++g)
    {
      const B_T::value_type &b = boundary(surface, g, boundary.OUT);
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
  }

  //-------------------------------------------------------------------------//
  // LEAKAGE RESPONSE
  //-------------------------------------------------------------------------//

  for (size_t surface = 0; surface < 4; ++surface)
  {
    size_t dim  = surface / 2;
    size_t dim0 = d_spatial_dim[dim][0];
    response->leakage_response(surface, index_i.nodal) = 0.0;
    for (size_t g = 0; g < d_material->number_groups(); ++g)
    {
      const B_T::value_type &bo = boundary(surface, g, boundary.OUT);
      const B_T::value_type &bi = boundary(surface, g, boundary.IN);

      for (size_t i = 0; i < bo.size(); ++i)
      {
        double dx = d_mesh->width(dim0, i);
        response->leakage_response(surface, index_i.nodal) +=
          dx * (B_V::value(bo, i) - B_V::value(bi, i));
      }
    }
  }

  //-------------------------------------------------------------------------//
  // FLUX RESPONSES
  //-------------------------------------------------------------------------//

  expand_flux(response, index_i);

}

//---------------------------------------------------------------------------//
template <>
void ResponseSourceDetran<detran::BoundaryDiffusion<detran::_3D> >::
expand(const Boundary_T    &boundary,
       SP_response          response,
       const ResponseIndex &index_i)
{
  using namespace detran;
  typedef BoundaryTraits<_3D> B_T;
  typedef BoundaryValue<_3D> B_V;
  //std::cout << " EXPANDING " << std::endl;

  //-------------------------------------------------------------------------//
  // CURRENT RESPONSE
  //-------------------------------------------------------------------------//

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
      const B_T::value_type &b = boundary(surface, g, boundary.OUT);

      // Expand the first spatial variable
      vec2_dbl R_s0m_s1(o_s0 + 1, vec_dbl(b.size(), 0.0));
      vec_dbl s0m(o_s0 + 1, 0);
      for (size_t i = 0; i < b.size(); ++i) // For all s1
      {
//        std::cout << " b.size=" << b[i].size()
//                  << " basis[" << surface << "][0].size "
//                  << d_basis_s[surface][0]->basis()->number_rows() << " "
//                  << d_basis_s[surface][0]->basis()->number_columns() << " "
//                  << " s0m size = " << s0m.size() << std::endl;
        d_basis_s[surface][0]->transform(b[i], s0m);
        for (size_t j = 0; j < s0m.size(); ++j)
          R_s0m_s1[j][i] = s0m[j];
        //THROW("lala");
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
    double sign = 1;
    if (surface == d_mesh->EAST  or
        surface == d_mesh->SOUTH or
        surface == d_mesh->BOTTOM)
    {
      sign = -1;
    }
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

  //-------------------------------------------------------------------------//
  // LEAKAGE RESPONSE
  //-------------------------------------------------------------------------//

  for (size_t surface = 0; surface < 6; ++surface)
  {
    size_t dim  = surface / 2;
    size_t dim0 = d_spatial_dim[dim][0];
    size_t dim1 = d_spatial_dim[dim][1];
    response->leakage_response(surface, index_i.nodal) = 0.0;
    for (size_t g = 0; g < d_material->number_groups(); ++g)
    {
      const B_T::value_type &bo = boundary(surface, g, boundary.OUT);
      const B_T::value_type &bi = boundary(surface, g, boundary.IN);

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

  //-------------------------------------------------------------------------//
  // FLUX RESPONSES
  //-------------------------------------------------------------------------//

  expand_flux(response, index_i);

}


//---------------------------------------------------------------------------//
// SN SPECIALIZATIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
template <>
void ResponseSourceDetran<detran::BoundarySN<detran::_1D> >::
set_boundary(Boundary_T &boundary, const ResponseIndex &index_i)
{
  //std::cout << "SETTING SOURCE: " << index_i << std::endl;
  using namespace detran;
  size_t octant = d_quadrature->incident_octant(index_i.surface)[0];
  for (size_t g = 0; g < d_material->number_groups(); ++g)
  {
    double P_e = (*d_basis_e[index_i.surface])(index_i.energy, g);
    for (size_t p = 0; p < d_quadrature->number_angles_octant(); ++p)
    {
      BoundaryTraits<_1D>::value_type
        &b = boundary(index_i.surface, octant, p, g);
      // \todo This is hard-coding an angular flux expansion
      double P_p = (*d_basis_p[index_i.surface])(index_i.polar, p);
      std::cout << " INCIDENT b[" << index_i.surface
                << " , " << p << "] = " << P_e * P_p << std::endl;
      BoundaryValue<_1D>::value(b) =  P_e * P_p;
    }
  }
}

//---------------------------------------------------------------------------//
template <>
void ResponseSourceDetran<detran::BoundarySN<detran::_1D> >::
expand(const Boundary_T    &boundary,
       SP_response          response,
       const ResponseIndex &index_i)
{

  //-------------------------------------------------------------------------//
  // CURRENT RESPONSE
  //-------------------------------------------------------------------------//
  //d_solver->state()->display();
  for (size_t surface = 0; surface < 2; ++surface)
  {
    size_t o_p = d_basis_p[surface]->order();
    size_t o_e = d_basis_e[surface]->order();
    size_t n_g = d_material->number_groups();
    size_t octant = d_quadrature->outgoing_octant(surface)[0];

    // Temporary response container
    vec2_dbl R(o_p + 1, vec_dbl(n_g, 0.0));
    vec_dbl Rp(o_p + 1, 0);

    // First expand in angle, [angle moments][energy groups]
    for (size_t g = 0; g < n_g; ++g)
    {
      vec_dbl psi_g(d_quadrature->number_angles_octant(), 0.0);
      for (size_t p = 0; p < d_quadrature->number_angles_octant(); ++p)
      {
        // \todo This hardcodes an expansion of the angular flux
        psi_g[p] = boundary(surface, octant, p, g);
        std::cout << " OUTGOING b[" << surface
                  << " , " << p << "] = " << psi_g[p] << std::endl;
      }
      d_basis_p[surface]->transform(psi_g, Rp);
      vec_dbl tmp(d_quadrature->number_angles_octant(), 0);
      d_basis_p[surface]->inverse(Rp, tmp);

      for (size_t p = 0; p < Rp.size(); ++p)
      {
        std::cout << " OUTGOING 2 b[" << surface
                  << " , " << p << "] = " << tmp[p] << std::endl;
        std::cout << " COEFvb[" << surface
                  << " , " << p << "] = " << Rp[p] << std::endl;
        R[p][g] = Rp[p];
      }
    }
    // Expand the result in energy, [angle moments][energy moments]
    for (size_t p = 0; p <= o_p; ++p)
    {
      vec_dbl R_g(o_e + 1, 0.0);
      d_basis_e[surface]->transform(R[p], R_g);
      R[p] = R_g;
    }
    // Fill the response container with only the *needed* values
    size_t nm = d_indexer->number_surface_moments(index_i.node, surface);
    for (size_t m = 0; m < nm; ++m)
    {
      ResponseIndex index_o = d_indexer->response_index(index_i.node, surface, m);
      response->boundary_response(index_o.nodal, index_i.nodal) =
        R[index_o.polar][index_o.energy];
    }
  }

  //-------------------------------------------------------------------------//
  // LEAKAGE RESPONSE
  //-------------------------------------------------------------------------//

  for (size_t surface = 0; surface < 2; ++surface)
  {
    size_t octant = d_quadrature->outgoing_octant(surface)[0];
    response->leakage_response(surface, index_i.nodal) = 0.0;
    for (size_t g = 0; g < d_material->number_groups(); ++g)
    {
      for (size_t p = 0; p < d_quadrature->number_angles_octant(); ++p)
      {
        const BoundaryTraits_T::value_type &psi_p =
          boundary(surface, octant, p, g);
        double w = d_quadrature->weight(p);
        double mu = d_quadrature->mu(0, p);
        response->leakage_response(surface, index_i.nodal) += mu * w * psi_p;
      }
    }
  }

  //-------------------------------------------------------------------------//
  // FLUX RESPONSES
  //-------------------------------------------------------------------------//

  expand_flux(response, index_i);

}

//---------------------------------------------------------------------------//
template <>
void ResponseSourceDetran<detran::BoundarySN<detran::_2D> >::
set_boundary(Boundary_T &boundary, const ResponseIndex &index)
{

}


//---------------------------------------------------------------------------//
template <>
void ResponseSourceDetran<detran::BoundarySN<detran::_2D> >::
expand(const Boundary_T    &boundary,
       SP_response          response,
       const ResponseIndex &index_i)
{

}


////---------------------------------------------------------------------------//
//template <>
//template <>
//void ResponseSourceDetran<detran::_3D>::
//set_boundary(detran::BoundarySN<detran::_3D>& boundary,
//             const ResponseIndex &index)
//{
//
//}
//
////---------------------------------------------------------------------------//
//template <>
//template <>
//void ResponseSourceDetran<detran::_3D>::
//expand(const detran::BoundarySN<detran::_3D>  &boundary,
//       SP_response                             response,
//       const ResponseIndex                     &index_i)
//{
//
//}
//
////---------------------------------------------------------------------------//
//// MOC SPECIALIZATIONS
////---------------------------------------------------------------------------//
//
////---------------------------------------------------------------------------//
//template <>
//template <>
//void ResponseSourceDetran<detran::_2D>::
//set_boundary(detran::BoundaryMOC<detran::_2D>& boundary,
//             const ResponseIndex &index)
//{
//  THROW("NOT IMPLEMENTED");
//}




} // end namespace erme_response

#endif // erme_response_RESPONSESOURCEDETRAN_T_HH_

//---------------------------------------------------------------------------//
//              end of file ResponseSourceDetran.t.hh
//---------------------------------------------------------------------------//
