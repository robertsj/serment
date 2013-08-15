//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ResponseSourceDetran.t.hh
 *  @brief ResponseSourceDetran
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef erme_response_RESPONSESOURCEDETRAN_T_HH_
#define erme_response_RESPONSESOURCEDETRAN_T_HH_

#include "ResponseSourceDetran.hh"
#include "boundary/BoundarySN.hh"
#include "boundary/BoundaryMOC.hh"
#include "boundary/BoundaryTraits.hh"
#include "boundary/FixedBoundary.hh"

#define DBOUT(c) std::cout << c << std::endl;

namespace erme_response
{

// \todo Need to consolidate code where possible

//----------------------------------------------------------------------------//
template <class B>
void ResponseSourceDetran<B>::set_boundary(const ResponseIndex &index)
{
  THROW("NOT IMPLEMENTED");
}
//----------------------------------------------------------------------------//
template <class B>
void ResponseSourceDetran<B>::expand_boundary(SP_response          response,
                                              const ResponseIndex &index_i)
{
  THROW("NOT IMPLEMENTED");
}

//----------------------------------------------------------------------------//
// SN SPECIALIZATIONS
//----------------------------------------------------------------------------//

#define FUCK(c) std::cout << c << std::endl;

//----------------------------------------------------------------------------//
// fixed boundary
template <>
void ResponseSourceDetran<detran::BoundarySN<detran::_1D> >::
set_boundary(const ResponseIndex &index_i)
{
  Boundary_T &B = *d_B;
  typedef detran::FixedBoundary<detran::_1D> Fixed_T;
  Fixed_T &BC = *dynamic_cast<Fixed_T*>(B.bc(index_i.surface).bp());
  size_t octant = d_quadrature->incident_octant(index_i.surface)[0];
  for (size_t g = 0; g < d_material->number_groups(); ++g)
  {
    double P_e = (*d_basis_e[index_i.surface])(g, index_i.energy);
    for (size_t p = 0; p < d_quadrature->number_angles_octant(); ++p)
    {


      double P_p = (*d_basis_p[index_i.surface])(p, index_i.polar);

      double val = P_e * P_p;

      if (!d_angular_flux) val /= d_quadrature->mu(0, p);
      BC(0, p, g) =  val;

      //B(index_i.surface, octant, p, g) = val;
    }
  }
}
// old boundary
//template <>
//void ResponseSourceDetran<detran::BoundarySN<detran::_1D> >::
//set_boundary(const ResponseIndex &index_i)
//{
//  Boundary_T &B = *d_B;
//  typedef detran::FixedBoundary<detran::_1D> Fixed_T;
//  size_t octant = d_quadrature->incident_octant(index_i.surface)[0];
//  for (size_t g = 0; g < d_material->number_groups(); ++g)
//  {
//    double P_e = (*d_basis_e[index_i.surface])(index_i.energy, g);
//    for (size_t p = 0; p < d_quadrature->number_angles_octant(); ++p)
//    {
//      double P_p = (*d_basis_p[index_i.surface])(index_i.polar, p);
//      double val = P_e * P_p;
//      if (!d_angular_flux) val /= d_quadrature->mu(0, p);
//      B(index_i.surface, octant, p, g) = val;
//    }
//  }
//  //B.display(true);
//}
//----------------------------------------------------------------------------//
template <>
void ResponseSourceDetran<detran::BoundarySN<detran::_1D> >::
expand_boundary(SP_response          response,
                const ResponseIndex &index_i)
{
  const Boundary_T &B = *d_B;

  for (size_t surface = 0; surface < 2; ++surface)
  {

    //------------------------------------------------------------------------//
    // BOUNDARY RESPONSE
    //------------------------------------------------------------------------//

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
        psi_g[p] = B(surface, octant, p, g);
        if (!d_angular_flux) psi_g[p] *= d_quadrature->mu(0, p);
      }
      d_basis_p[surface]->transform(psi_g, Rp);
      for (size_t p = 0; p < Rp.size(); ++p)
      {
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

    //------------------------------------------------------------------------//
    // LEAKAGE
    //------------------------------------------------------------------------//

    response->leakage_response(surface, index_i.nodal) = 0.0;
    for (size_t g = 0; g < d_material->number_groups(); ++g)
    {
      for (size_t p = 0; p < d_quadrature->number_angles_octant(); ++p)
      {
        double psi_p = B(surface, octant, p, g);
        double w = d_quadrature->weight(p);
        double mu = d_quadrature->mu(0, p);
        response->leakage_response(surface, index_i.nodal) += mu * w * psi_p;
      }
    }
  }
//  response->display();
//  THROW("lala");
}

//----------------------------------------------------------------------------//
template <>
void ResponseSourceDetran<detran::BoundarySN<detran::_2D> >::
set_boundary(const ResponseIndex &index_i)
{
  Boundary_T &B = *d_B;
  typedef detran::FixedBoundary<detran::_2D> Fixed_T;
  Fixed_T &BC = *dynamic_cast<Fixed_T*>(B.bc(index_i.surface).bp());

  double sign = 1.0;
  if ((index_i.surface == d_mesh->WEST   or
       index_i.surface == d_mesh->NORTH) and
      (index_i.space0) % 2 )
  {
    sign = -1.0;
  }

  SP_productquadrature q = d_quadrature;

  size_t dim  = index_i.surface / 2;
  size_t dim0 = d_spatial_dim[dim][0];

  int na = q->number_azimuths_octant();

  //                    octant 0        octant 1
  int az[4][2][3] = { {{na-1, -1, -1}, {0,     1, na}},    // west
                      {{na-1, -1, -1}, {0,     1, na}},    // east
                      {{0,     1, na}, {na-1, -1, -1}},    // south
                      {{0,     1, na}, {na-1, -1, -1}} };  // north

  /*
   *  Need a clean way to loop over azimuths in right order.  Detran
   *  product quads are guaranteed to give angles such that for each
   *  polar level, mu goes from 0 : 1.  In that way, octant 1 is
   *  correctly ordered for incident on the left (going left to right),
   *  but octant 4 must be reversed in azimuth.
   */
  for (size_t g = 0; g < d_material->number_groups(); ++g)
  {
    //DBOUT(" g = " << g)
    double P_e = (*d_basis_e[index_i.surface])(index_i.energy, g);
    size_t a = 0;
    for (size_t oo = 0; oo < 2; ++oo)
    {
      size_t o = q->incident_octant(index_i.surface)[oo];
      //DBOUT("   o = " << o)
      int a0 = az[index_i.surface][oo][1];
      for (size_t aa = 0; aa < q->number_azimuths_octant(); ++aa, ++a)
      {
        size_t aaa = aa;
        if (a0 == -1) aaa = q->number_azimuths_octant() - aa - 1;
        //DBOUT("     aaa = " << aaa << " " << q->number_azimuths_octant() << " " << a0)

        double P_a = (*d_basis_a[index_i.surface])(index_i.azimuth, a);
        for (size_t p = 0; p < q->number_polar_octant(); ++p)
        {
          //DBOUT("       p = " << p)
          double P_p = (*d_basis_p[index_i.surface])(index_i.polar, p);
          size_t angle = q->angle(aaa, p);
          BoundaryTraits_T::value_type &b = BC(oo, angle, g);
          for (size_t i = 0; i < d_mesh->number_cells(dim0); ++i)
          {
            double P_s0 = (*d_basis_s[index_i.surface][0])(index_i.space0, i);
            double val = sign * P_s0 * P_p * P_a * P_e;
            if (!d_angular_flux) val /= q->cosines(dim)[angle];
            BoundaryValue_T::value(b, i) = val;
          }
        }
      }
    }
  }
}

// OLD
//template <>
//void ResponseSourceDetran<detran::BoundarySN<detran::_2D> >::
//set_boundary(const ResponseIndex &index_i)
//{
//  Boundary_T &B = *d_B;
//  double sign = 1.0;
//  if ((index_i.surface == d_mesh->WEST   or
//       index_i.surface == d_mesh->NORTH) and
//      (index_i.space0) % 2 )
//  {
//    sign = -1.0;
//  }
//
//  SP_productquadrature q = d_quadrature;
//
//  size_t dim  = index_i.surface / 2;
//  size_t dim0 = d_spatial_dim[dim][0];
//
//  int na = q->number_azimuths_octant();
//
//  //                    octant 0        octant 1
//  int az[4][2][3] = { {{na-1, -1, -1}, {0,     1, na}},    // west
//                      {{na-1, -1, -1}, {0,     1, na}},    // east
//                      {{0,     1, na}, {na-1, -1, -1}},    // south
//                      {{0,     1, na}, {na-1, -1, -1}} };  // north
//
//  /*
//   *  Need a clean way to loop over azimuths in right order.  Detran
//   *  product quads are guaranteed to give angles such that for each
//   *  polar level, mu goes from 0 : 1.  In that way, octant 1 is
//   *  correctly ordered for incident on the left (going left to right),
//   *  but octant 4 must be reversed in azimuth.
//   */
//  for (size_t g = 0; g < d_material->number_groups(); ++g)
//  {
//    //std::cout << " g = " << g << std::endl;
//    double P_e = (*d_basis_e[index_i.surface])(index_i.energy, g);
//    size_t a = 0;
//    for (size_t oo = 0; oo < 2; ++oo)
//    {
//      size_t o = q->incident_octant(index_i.surface)[oo];
//      //std::cout << "   o = " << o << std::endl;
//      int a0 = az[index_i.surface][oo][1];
//      for (size_t aa = 0; aa < q->number_azimuths_octant(); ++aa, ++a)
//      {
//        size_t aaa = aa;
//        if (a0 == -1) aaa = q->number_azimuths_octant() - aa - 1;
//
//        //std::cout << "     aaa = " << aaa << " " << q->number_azimuths_octant() << " " << a0 << std::endl;
//        double P_a = (*d_basis_a[index_i.surface])(index_i.azimuth, a);
//        for (size_t p = 0; p < q->number_polar_octant(); ++p)
//        {
//          //std::cout << "       p = " << p << std::endl;
//          double P_p = (*d_basis_p[index_i.surface])(index_i.polar, p);
//          size_t angle = q->angle(aaa, p);
//          BoundaryTraits_T::value_type &b = B(index_i.surface, o, angle, g);
//          for (size_t i = 0; i < d_mesh->number_cells(dim0); ++i)
//          {
//            double P_s0 = (*d_basis_s[index_i.surface][0])(index_i.space0, i);
//            double val = sign * P_s0 * P_p * P_a * P_e;
//            if (!d_angular_flux) val /= q->cosines(dim)[angle];
//            BoundaryValue_T::value(b, i) = val;
//          }
//        }
//      }
//    }
//  }
//}

//----------------------------------------------------------------------------//
template <>
void ResponseSourceDetran<detran::BoundarySN<detran::_2D> >::
expand_boundary(SP_response          response,
                const ResponseIndex &index_i)
{

  const Boundary_T &B = *d_B;
  SP_productquadrature q = d_quadrature;

  int na = q->number_azimuths_octant();

  //                    octant 0        octant 1
  int az[4][2][3] = { {{na-1, -1, -1}, {0,     1, na}},    // west
                      {{na-1, -1, -1}, {0,     1, na}},    // east
                      {{0,     1, na}, {na-1, -1, -1}},    // south
                      {{0,     1, na}, {na-1, -1, -1}} };  // north

  for (size_t surface = 0; surface < 4; ++surface)
  {
    size_t dim  = surface / 2;
    size_t dim0 = d_spatial_dim[dim][0];

    double sign = 1;
    if (surface == d_mesh->EAST or surface == d_mesh->SOUTH)
    {
      sign = -1;
    }

    //------------------------------------------------------------------------//
    // BOUNDARY RESPONSE
    //------------------------------------------------------------------------//

    size_t o_e = d_basis_e[surface]->order();
    size_t o_a = d_basis_a[surface]->order();
    size_t o_p = d_basis_p[surface]->order();
    size_t o_s = d_basis_s[surface][0]->order();
    size_t n_g = d_material->number_groups();

    // Temporary response containers
    vec4_dbl R(o_s + 1,
               vec3_dbl(o_p + 1,
                        vec2_dbl(o_a + 1,
                                 vec_dbl(n_g, 0.0))));

    // First expand in angle, [angle moments][energy groups]
    for (size_t g = 0; g < n_g; ++g)
    {

      vec3_dbl R_s_p_a(o_s + 1,
                       vec2_dbl(o_p + 1,
                                vec_dbl(2 * q->number_azimuths_octant(), 0)));

      // Azimuth block
      {
        size_t a = 0;
        for (size_t oo = 0; oo < 2; ++oo)
        {
          size_t o = q->outgoing_octant(surface)[oo];
          int a0 = az[surface][oo][1];

          for (size_t aa = 0; aa < q->number_azimuths_octant(); ++aa, ++a)
          {
            size_t aaa = aa;
            if (a0 == -1) aaa = q->number_azimuths_octant() - aa - 1;

            vec2_dbl R_s_p(o_s + 1, vec_dbl(q->number_polar_octant(), 0));

            for (size_t p = 0; p < q->number_polar_octant(); ++p)
            {
              size_t angle = q->angle(aaa, p);

              const BoundaryTraits_T::value_type &b = B(surface, o, angle, g);
              vec_dbl R_s(o_s + 1, 0);
              //std::cout << " surface=" << surface << " b=" << b[0] << std::endl;
              // EXPAND S
              d_basis_s[surface][0]->transform(b, R_s);
              for (size_t s = 0; s < R_s.size(); ++s)
              {
                if (!d_angular_flux) R_s[s] *= q->cosines(dim)[angle];
                R_s_p[s][p] = R_s[s];
              }

            } // end p

            // EXPAND P
            vec_dbl R_p(o_p + 1, 0);
            for (size_t s = 0; s < R_s_p.size(); ++s)
            {
              d_basis_p[surface]->transform(R_s_p[s], R_p);
              for (size_t p = 0; p < R_p.size(); ++p)
                R_s_p_a[s][p][a] = R_p[p];
            }

          } // end a
        } // end o
      } // end azimuth

      // EXPAND A
      vec_dbl R_a(o_a + 1, 0);
      for (size_t s = 0; s < R_s_p_a.size(); ++s)
      {
        for (size_t p = 0; p < R_s_p_a[s].size(); ++p)
        {
          d_basis_a[surface]->transform(R_s_p_a[s][p], R_a);
          for (size_t a = 0; a < R_a.size(); ++a)
            R[s][p][a][g] = R_a[a];
        }
      }

    } // end g

    // EXPAND G
    vec_dbl R_g(o_e + 1, 0);
    for (size_t s = 0; s < R.size(); ++s)
    {
      for (size_t p = 0; p < R[s].size(); ++p)
      {
        for (size_t a = 0; a < R[s][p].size(); ++a)
        {
          d_basis_e[surface]->transform(R[s][p][a], R_g);
          R[s][p][a] = R_g;
        }
      }
    }

    // Fill the response container with only the *needed* values
    size_t nm = d_indexer->number_surface_moments(index_i.node, surface);
    for (size_t m = 0; m < nm; ++m)
    {
      ResponseIndex index_o =
        d_indexer->response_index(index_i.node, surface, m);
      double coef = std::pow(sign, index_o.space0);
      response->boundary_response(index_o.nodal, index_i.nodal) =
        coef * R[index_o.space0][index_o.polar][index_o.azimuth][index_o.energy];
    }

    //------------------------------------------------------------------------//
    // LEAKAGE
    //------------------------------------------------------------------------//

    response->leakage_response(surface, index_i.nodal) = 0.0;
    for (size_t g = 0; g < d_material->number_groups(); ++g)
    {
      for (size_t oo = 0; oo < 2; ++oo)
      {
        size_t o = d_quadrature->outgoing_octant(surface)[oo];
        for (size_t a = 0; a < d_quadrature->number_angles_octant(); ++a)
        {
          const BoundaryTraits_T::value_type &psi = B(surface, o, a, g);
          double w = d_quadrature->weight(a);
          double mu = d_quadrature->cosines(dim)[a];
          for (size_t s = 0; s < d_mesh->number_cells(dim0); ++s)
          {
            double dx = d_mesh->width(dim0, s);
            response->leakage_response(surface, index_i.nodal) +=
              mu * w * dx * BoundaryValue_T::value(psi, s);
          }
        }
      }
    }

  } // end surface

}


////----------------------------------------------------------------------------//
//template <>
//template <>
//void ResponseSourceDetran<detran::_3D>::
//set_boundary(detran::BoundarySN<detran::_3D>& boundary,
//             const ResponseIndex &index)
//{
//
//}
//
////----------------------------------------------------------------------------//
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
////----------------------------------------------------------------------------//
//// MOC SPECIALIZATIONS
////----------------------------------------------------------------------------//
//
////----------------------------------------------------------------------------//
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

//----------------------------------------------------------------------------//
//              end of file ResponseSourceDetran.t.hh
//----------------------------------------------------------------------------//
