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

//---------------------------------------------------------------------------//
template <class D>
template <class B>
void ResponseSourceDetran<D>::set_boundary(B& boundary,
                                           const ResponseIndex &index)
{
  THROW("NOT IMPLEMENTED");
}

//---------------------------------------------------------------------------//
// DIFFUSION SPECIALIZATIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
template <>
template <>
void ResponseSourceDetran<detran::_1D>::
set_boundary(detran::BoundaryDiffusion<detran::_1D>& boundary,
             const ResponseIndex &index)
{
  using namespace detran;
  for (size_t g = 0; g < d_material->number_groups(); ++g)
  {
    std::cout << " g = " << g << std::endl;
    BoundaryTraits<_1D>::value_type
      &b = boundary(index.surface, g, boundary.IN);
    Assert(d_basis_e[index.surface]);
    BoundaryValue<_1D>::value(b, 0, 0) =
      (*d_basis_e[index.surface])(index.energy, g);
  }
}

//---------------------------------------------------------------------------//
template <>
template <>
void ResponseSourceDetran<detran::_2D>::
set_boundary(detran::BoundaryDiffusion<detran::_2D>& boundary,
             const ResponseIndex &index)
{
  using namespace detran;
  double sign = 1.0;
  if ((index.surface == d_mesh->WEST or
       index.surface == d_mesh->NORTH) and index.space0 % 2)
  {
    sign = -1.0;
  }
  for (size_t g = 0; g < d_material->number_groups(); ++g)
  {
    double P_e = sign * (*d_basis_e[index.surface])(index.energy, g);
    size_t dim  = index.surface / 2;
    size_t dim0 = d_spatial_dim[dim][0];
    BoundaryTraits<_2D>::value_type
      &b = boundary(index.surface, g, boundary.IN);
    for (size_t i = 0; i < d_mesh->number_cells(dim0); ++i)
    {
      double P_s0 = (*d_basis_s[index.surface][0])(index.space0, i);
      BoundaryValue<_2D>::value(b, i, 0) =  P_e * P_s0;
    }
  }
}

//---------------------------------------------------------------------------//
template <>
template <>
void ResponseSourceDetran<detran::_3D>::
set_boundary(detran::BoundaryDiffusion<detran::_3D>& boundary,
             const ResponseIndex &index)
{
  using namespace detran;
  for (size_t g = 0; g < d_material->number_groups(); ++g)
  {
    double P_e = (*d_basis_e[index.surface])(index.energy, g);
    size_t dim  = index.surface / 2;
    size_t dim0 = d_spatial_dim[dim][0];
    size_t dim1 = d_spatial_dim[dim][1];
    BoundaryTraits<_3D>::value_type
      &b = boundary(g, index.surface, boundary.IN);
    for (size_t i = 0; i < d_mesh->number_cells(dim0); ++i)
    {
      double P_s0 = (*d_basis_s[index.surface][0])(index.space0, i);
      for (size_t j = 0; j < d_mesh->number_cells(dim1); ++j)
      {
        double P_s1 = (*d_basis_s[index.surface][1])(index.space1, j);
        BoundaryValue<_3D>::value(b, i, j) = P_e * P_s0 * P_s1;
      }
    }
  }
}

//---------------------------------------------------------------------------//
template <>
template <>
void ResponseSourceDetran<detran::_1D>::
expand(const detran::BoundaryDiffusion<detran::_1D>  &boundary,
       SP_response                                    response,
       const ResponseIndex                           &index_i)
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
    for (size_t g = 0; g < d_material->number_groups(); ++g)
    {
      const B_T::value_type &bo = boundary(g, surface, boundary.OUT);
      const B_T::value_type &bi = boundary(g, surface, boundary.IN);
      response->leakage_response(surface, index_i.nodal) +=
        (B_V::value(bo, 0, 0) - B_V::value(bi, 0, 0));
    }
  }

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
template <>
template <>
void ResponseSourceDetran<detran::_2D>::
expand(const detran::BoundaryDiffusion<detran::_2D>  &boundary,
       SP_response                                    response,
       const ResponseIndex                           &index_i)
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
        R[s][g] = Rs[s];
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
//      if (index_o.nodal == 1 and index_i.nodal == 1)
//      {
//        std::cout << index_i << index_o << " poo = " << poo << std::endl;
//        THROW("lalala");
//      }
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
template <>
template <>
void ResponseSourceDetran<detran::_3D>::
expand(const detran::BoundaryDiffusion<detran::_3D>  &boundary,
       SP_response                                    response,
       const ResponseIndex                           &index_i)
{
  using namespace detran;
  typedef BoundaryTraits<_3D> B_T;
  typedef BoundaryValue<_3D> B_V;

  for (size_t surface = 0; surface < 6; ++surface)
  {

  }
}


//---------------------------------------------------------------------------//
// SN SPECIALIZATIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
template <>
template <>
void ResponseSourceDetran<detran::_1D>::
set_boundary(detran::BoundarySN<detran::_1D>& boundary,
             const ResponseIndex &index)
{

}

//---------------------------------------------------------------------------//
template <>
template <>
void ResponseSourceDetran<detran::_2D>::
set_boundary(detran::BoundarySN<detran::_2D>& boundary,
             const ResponseIndex &index)
{

}

//---------------------------------------------------------------------------//
template <>
template <>
void ResponseSourceDetran<detran::_3D>::
set_boundary(detran::BoundarySN<detran::_3D>& boundary,
             const ResponseIndex &index)
{

}

//---------------------------------------------------------------------------//
// MOC SPECIALIZATIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
template <>
template <>
void ResponseSourceDetran<detran::_2D>::
set_boundary(detran::BoundaryMOC<detran::_2D>& boundary,
             const ResponseIndex &index)
{
  THROW("NOT IMPLEMENTED");
}



} // end namespace erme_response

#endif // erme_response_RESPONSESOURCEDETRAN_T_HH_

//---------------------------------------------------------------------------//
//              end of file ResponseSourceDetran.t.hh
//---------------------------------------------------------------------------//
