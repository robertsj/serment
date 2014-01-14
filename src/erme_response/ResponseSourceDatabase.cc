//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  ResponseSourceDatabase.cc
 *  @brief ResponseSourceDatabase member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#include "ResponseSourceDatabase.hh"

namespace erme_response
{

//---------------------------------------------------------------------------//
ResponseSourceDatabase::
ResponseSourceDatabase(SP_node node, SP_indexer indexer)
  : ResponseSource(node, indexer)
  , d_rfdb(ResponseDatabase::Instance())
{
  /* ... */
}

//---------------------------------------------------------------------------//
ResponseSourceDatabase::~ResponseSourceDatabase()
{
  /* ... */
}

//---------------------------------------------------------------------------//
void ResponseSourceDatabase::
compute(SP_response response, const ResponseIndex &index)
{
  d_rfdb->get(d_node->name(), response, index, d_keff);
}

} // end namespace erme_response

//----------------------------------------------------------------------------//
//              end of file ResponseSourceDatabase.cc
//----------------------------------------------------------------------------//
