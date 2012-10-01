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

namespace erme_response
{

inline ResponseSourceFactory::SP_source
ResponseSourceFactory::build_detran(SP_node node)
{
  SP_source s;
  return s;
}

} // end namespace erme_response

#endif /* RESPONSESOURCEFACTORYDETRAN_HH_ */
