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
                                           ResponseIndex index)
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
             ResponseIndex index)
{
  using namespace detran;
  for (size_t g = 0; g < d_material->number_groups(); ++g)
  {
    BoundaryTraits<_1D>::value_type
      &b = boundary(g, index.surface, boundary.IN);
    BoundaryValue<_1D>::value(b, 0, 0) =
      (*d_basis_e[index.surface]->basis())(index.energy, g);
  }
}

//---------------------------------------------------------------------------//
template <>
template <>
void ResponseSourceDetran<detran::_2D>::
set_boundary(detran::BoundaryDiffusion<detran::_2D>& boundary,
             ResponseIndex index)
{
  using namespace detran;
}

//---------------------------------------------------------------------------//
template <>
template <>
void ResponseSourceDetran<detran::_3D>::
set_boundary(detran::BoundaryDiffusion<detran::_3D>& boundary,
             ResponseIndex index)
{

}

//---------------------------------------------------------------------------//
// SN SPECIALIZATIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
template <>
template <>
void ResponseSourceDetran<detran::_1D>::
set_boundary(detran::BoundarySN<detran::_1D>& boundary,
             ResponseIndex index)
{

}

//---------------------------------------------------------------------------//
template <>
template <>
void ResponseSourceDetran<detran::_2D>::
set_boundary(detran::BoundarySN<detran::_2D>& boundary,
             ResponseIndex index)
{

}

//---------------------------------------------------------------------------//
template <>
template <>
void ResponseSourceDetran<detran::_3D>::
set_boundary(detran::BoundarySN<detran::_3D>& boundary,
             ResponseIndex index)
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
             ResponseIndex index)
{

}



} // end namespace erme_response

#endif // erme_response_RESPONSESOURCEDETRAN_T_HH_

//---------------------------------------------------------------------------//
//              end of file ResponseSourceDetran.t.hh
//---------------------------------------------------------------------------//
