//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ResponseSourceFactoryDatabase.hh
 *  @author robertsj
 *  @date   Oct 1, 2012
 *  @brief  ResponseSourceFactoryDatabase class definition.
 */
//---------------------------------------------------------------------------//

#ifndef erme_response_RESPONSESOURCEFACTORYDATABASE_HH_
#define erme_response_RESPONSESOURCEFACTORYDATABASE_HH_

#include "erme_geometry/CartesianNode.hh"

namespace erme_response
{

//---------------------------------------------------------------------------//
inline ResponseSourceFactory::SP_source
ResponseSourceFactory::build_database(SP_node node, SP_indexer indexer)
{
  SP_source s(new ResponseSourceDatabase(node, indexer));
  return s;
}

} // end namespace erme_response

#endif /* erme_response_RESPONSESOURCEFACTORYDATABASE_HH_ */
