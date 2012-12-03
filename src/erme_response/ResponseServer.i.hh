//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ResponseServer.i.hh
 *  @author robertsj
 *  @date   Aug 31, 2012
 *  @brief  ResponseServer inline member definitions
 */
//---------------------------------------------------------------------------//

#ifndef erme_response_RESPONSESERVER_I_HH_
#define erme_response_RESPONSESERVER_I_HH_

#include "comm/Comm.hh"

namespace erme_response
{

//---------------------------------------------------------------------------//
inline ResponseServer::SP_response ResponseServer::response(size_t node_l)
{
  // Preconditions \todo not an optimal sequence of indexing...
  size_t node_g  = d_nodes->global_index(node_l);
  size_t node_ug = d_nodes->unique_global_index(node_g);
  size_t node_ul = d_nodes->unique_local_index(node_ug);
  Require(node_ul < d_responses.size());

  return d_responses[node_ul];
}

} // end namespace erme_response

#endif /* erme_response_RESPONSESERVER_I_HH_ */
