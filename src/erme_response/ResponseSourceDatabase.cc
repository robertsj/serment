//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ResponseSourceDatabase.cc
 *  @author robertsj
 *  @date   Oct 1, 2012
 *  @brief  ResponseSourceDatabase class definition.
 */
//---------------------------------------------------------------------------//

#include "ResponseSourceDatabase.hh"

namespace erme_response
{

ResponseSourceDatabase::ResponseSourceDatabase(SP_node node)
  : ResponseSource(node)
  , d_rfdb(ResponseDatabase::Instance())
{
  /* ... */
}

ResponseSourceDatabase::~ResponseSourceDatabase()
{
  /* ... */
}

void ResponseSourceDatabase::
compute(SP_response response, ResponseIndex index)
{
  d_rfdb->get(d_node->name(), response, index, d_keff);
}

} // end namespace erme_response
