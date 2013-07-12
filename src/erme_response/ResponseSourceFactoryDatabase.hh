//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ResponseSourceFactoryDatabase.hh
 *  @brief ResponseSourceFactoryDatabase class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef erme_response_RESPONSESOURCEFACTORYDATABASE_HH_
#define erme_response_RESPONSESOURCEFACTORYDATABASE_HH_

#include "erme_geometry/CartesianNode.hh"

namespace erme_response
{

//----------------------------------------------------------------------------//
inline ResponseSourceFactory::SP_source
ResponseSourceFactory::build_database(SP_node node, SP_indexer indexer)
{
  SP_source s(new ResponseSourceDatabase(node, indexer));
  return s;
}

} // end namespace erme_response

#endif /* erme_response_RESPONSESOURCEFACTORYDATABASE_HH_ */

//----------------------------------------------------------------------------//
//              end of file ResponseSourceFactoryDatabase.hh
//----------------------------------------------------------------------------//
