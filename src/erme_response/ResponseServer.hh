//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseServer.hh
 * \brief  ResponseServer class definition
 * \author Jeremy Roberts
 * \date   Aug 28, 2012
 */
//---------------------------------------------------------------------------//

#ifndef RESPONSESERVER_HH_
#define RESPONSESERVER_HH_

#include "NodeResponse.hh"
#include "DBC.hh"
#include "SP.hh"
#include <vector>

namespace erme_response
{

/*!
 *  \class ResponseServer
 *  \brief Serve nodal responses to clients
 */
class ResponseServer
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran::SP<ResponseServer>  SP_server;
  typedef NodeResponse::SP_response   SP_response;
  typedef std::vector<SP_response>    vec_response;
  typedef unsigned int                size_t;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   */
  ResponseServer();

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Node response functions
  vec_response d_responses;


  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

};



} // end namespace erme_response

#endif // RESPONSESERVER_HH_ 

//---------------------------------------------------------------------------//
//              end of file ResponseServer.hh
//---------------------------------------------------------------------------//
