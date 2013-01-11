//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ResponseSourceFactoryDetran.hh
 *  @author robertsj
 *  @date   Aug 31, 2012
 *  @brief  ResponseSourceFactoryDetran build specialization
 */
//---------------------------------------------------------------------------//

#ifndef RESPONSESOURCEFACTORYDETRAN_HH_
#define RESPONSESOURCEFACTORYDETRAN_HH_

#include "erme_geometry/CartesianNodeDetran.hh"
#include "transport/DimensionTraits.hh"

namespace erme_response
{

//---------------------------------------------------------------------------//
inline ResponseSourceFactory::SP_source
ResponseSourceFactory::build_detran(SP_node node)
{
  erme_geometry::CartesianNodeDetran::SP_node detran_node = node;
  Require(detran_node->mesh());
  size_t d = detran_node->mesh()->dimension();

  SP_source s;
  if (d == 0)
    s = new ResponseSourceDetran<detran::_1D>(node);
  else if (d == 1)
    s = new ResponseSourceDetran<detran::_2D>(node);
  else
    s = new ResponseSourceDetran<detran::_3D>(node);
  return s;
}

} // end namespace erme_response

#endif /* RESPONSESOURCEFACTORYDETRAN_HH_ */
