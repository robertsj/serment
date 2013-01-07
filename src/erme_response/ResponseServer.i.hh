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
  // Preconditions
  Require(node_l < d_nodes->number_local_nodes());

  size_t node_g  = d_nodes->global_index_from_local(node_l);
  size_t node_ug = d_nodes->unique_global_index_from_global(node_g);
  size_t node_ul = d_nodes->unique_local_index_from_unique_global(node_ug);

//  std::cout << " node_l=" << node_l
//            << " node_g=" << node_g
//            << " node_ug=" << node_ug
//            << " node_ul=" << node_ul << std::endl;

  // Postconditions
  Ensure(node_ul < d_responses.size());
  return d_responses[node_ul];
}

} // end namespace erme_response

#endif /* erme_response_RESPONSESERVER_I_HH_ */
