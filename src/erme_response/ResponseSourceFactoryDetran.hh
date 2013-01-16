//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ResponseSourceFactoryDetran.hh
 *  @author robertsj
 *  @date   Aug 31, 2012
 *  @brief  ResponseSourceFactoryDetran build specialization
 */
//---------------------------------------------------------------------------//

#ifndef erme_response_RESPONSESOURCEFACTORYDETRAN_HH_
#define erme_response_RESPONSESOURCEFACTORYDETRAN_HH_

#include "erme_geometry/CartesianNodeDetran.hh"
#include "transport/DimensionTraits.hh"

namespace erme_response
{

//---------------------------------------------------------------------------//
inline ResponseSourceFactory::SP_source
ResponseSourceFactory::build_detran(SP_node node, SP_indexer indexer)
{
  erme_geometry::CartesianNodeDetran::SP_node detran_node = node;
  Require(detran_node->mesh());
  size_t d = detran_node->mesh()->dimension();

  SP_source s;
  if (d == 0)
    s = new ResponseSourceDetran<detran::_1D>(node, indexer);
  else if (d == 1)
    s = new ResponseSourceDetran<detran::_2D>(node, indexer);
  else
    s = new ResponseSourceDetran<detran::_3D>(node, indexer);
  return s;
}

} // end namespace erme_response

#endif /* erme_response_RESPONSESOURCEFACTORYDETRAN_HH_ */
