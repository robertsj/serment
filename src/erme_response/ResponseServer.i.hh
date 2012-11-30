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
inline ResponseServer::SP_response ResponseServer::response(size_t node)
{
  // Preconditions
  Require(node < d_responses.size());

  return d_responses[node];
}

} // end namespace erme_response

#endif /* erme_response_RESPONSESERVER_I_HH_ */
