//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ResponseSourceFactoryDummy.hh
 *  @author robertsj
 *  @date   Aug 31, 2012
 *  @brief  ResponseSourceFactory build specialization for Dummy nodes
 */
//---------------------------------------------------------------------------//

#ifndef erme_response_RESPONSESOURCEFACTORYDUMMY_HH_
#define erme_response_RESPONSESOURCEFACTORYDUMMY_HH_

#include "erme_geometry/DummyNode.hh"
#include "ResponseSourceDummy.hh"

namespace erme_response
{

//---------------------------------------------------------------------------//
inline ResponseSourceFactory::SP_source
ResponseSourceFactory::build_dummy(SP_node node, SP_indexer indexer)
{
  SP_source s(new ResponseSourceDummy(node, indexer));
  return s;
}

} // end namespace erme_response

#endif /* erme_response_RESPONSESOURCEFACTORYDUMMY_HH_ */
