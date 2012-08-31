//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseServer.cc
 * \author robertsj
 * \date   Aug 31, 2012
 * \brief  ResponseServer class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#include "ResponseServer.hh"
#include "ResponseSourceFactory.hh"

namespace erme_response
{

ResponseServer::ResponseServer(erme_geometry::NodeList &nodes, ResponseIndexer &indexer)
  : d_nodes(nodes)
  , d_indexer(indexer)
  , d_sources(nodes.number_local_nodes())
{

  ResponseSourceFactory builder;
  for (size_t n = 0; n < d_sources.size(); n++)
  {
    size_t n_global = nodes.global_index(n);
    d_sources[n] = builder.build(nodes.node(n_global));
    Assert(d_sources[n]);
  }

}

void ResponseServer::update(const double keff)
{

  // Broadcast keff

  // Update sources

  // Update responses

}

} // end namespace erme_response


