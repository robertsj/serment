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

#include "erme_geometry/NodeList.hh"
#include "ResponseIndexer.hh"
#include "NodeResponse.hh"
#include "ResponseSource.hh"
#include "DBC.hh"
#include "SP.hh"
#include <vector>

namespace erme_response
{

/*!
 *  \class ResponseServer
 *  \brief Serve nodal responses to clients
 *
 *  A ResponseServer lives on a local communicator.  A server is in charge
 *  of one or more Node objects.  The nodal responses are produced by a
 *  ResponseSource that solves the local problems.
 *
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
  typedef ResponseSource::SP_source   SP_source;
  typedef std::vector<SP_source>      vec_source;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   */
  ResponseServer(erme_geometry::NodeList &nodes, ResponseIndexer &indexer);

  /// Update the eigenvalue and compute the new responses
  void update(const double keff);

  /// Return a nodal response
  SP_response response(size_t node);

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Nodes
  erme_geometry::NodeList& d_nodes;

  /// Indexer
  ResponseIndexer& d_indexer;

  /// Sources
  vec_source d_sources;

  /// Node response functions
  vec_response d_responses;


  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Compute responses by dividing work explicitly among processes
   *
   *  Each response is treated the same, when in fact some responses (or
   *  nodes) might be inherently expensive.  A reduction is used to get
   *  the data, which probably is inefficient for lots of nodes and/or
   *  lots of response data.
   *
   */
  void update_explicit_work_share();

  /*!
   *  \brief Compute responses by using a self-scheduling master-slave
   *         approach
   */
  void update_master_slave()
  {
    THROW("not yet done");
  }

};

} // end namespace erme_response

//---------------------------------------------------------------------------//
// INLINE MEMBER DEFINITIONS
//---------------------------------------------------------------------------//

#include "ResponseServer.i.hh"

#endif // RESPONSESERVER_HH_ 

//---------------------------------------------------------------------------//
//              end of file ResponseServer.hh
//---------------------------------------------------------------------------//
