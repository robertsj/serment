//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseServer.i.hh
 * \author robertsj
 * \date   Aug 31, 2012
 * \brief  ResponseServer inline member definitions
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef RESPONSESERVER_I_HH_
#define RESPONSESERVER_I_HH_

#include "comm/Comm.hh"

namespace erme_response
{


inline ResponseServer::SP_response ResponseServer::response(size_t node)
{
  return d_responses[node];
}


} // end namespace erme_response

#endif /* RESPONSESERVER_I_HH_ */
