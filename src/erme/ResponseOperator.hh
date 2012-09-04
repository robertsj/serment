//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseOperator.hh
 * \brief  ResponseOperator class definition
 * \author Jeremy Roberts
 * \date   Aug 24, 2012
 */
//---------------------------------------------------------------------------//

#ifndef RESPONSEOPERATOR_HH_
#define RESPONSEOPERATOR_HH_

#include "erme_geometry/NodeList.hh"
#include "erme_response/ResponseIndexer.hh"
#include "erme_response/ResponseServer.hh"

#include "SP.hh"

namespace erme
{

/*!
 *  \class ResponseOperator
 *  \brief Base class for response operators
 */
class ResponseOperator
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef erme_geometry::NodeList::SP_nodelist        SP_nodelist;
  typedef erme_response::ResponseIndexer::SP_indexer  SP_indexer;
  typedef erme_response::ResponseServer::SP_server    SP_server;
  typedef erme_response::ResponseServer::SP_response  SP_response;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   *  \param indexer  Pointer to response indexer
   *  \param server   Pointer to response server
   */
  ResponseOperator(SP_nodelist nodes, SP_indexer indexer, SP_server server)
    : d_nodes(nodes)
    , d_indexer(indexer)
    , d_server(server)
  {
    // Preconditions
    Require(d_nodes);
    Require(d_indexer);
    Require(d_server);
  }

  /// Virtual Destructor
  virtual ~ResponseOperator(){}

  /*!
   *  \brief Update responses for new eigenvalue
   *
   *  Each operator must implement this for its
   *  own operator structure.
   *
   *  \param keff   New k-eigenvalue
   */
  virtual void update(const double keff) = 0;

protected:

  //-------------------------------------------------------------------------//
  // PROTECTED DATA
  //-------------------------------------------------------------------------//

  /// Node list
  SP_nodelist d_nodes;

  /// Response indexer
  SP_indexer d_indexer;

  /// Response server
  SP_server d_server;

};

} // end namespace erme

#endif // RESPONSEOPERATOR_HH_ 

//---------------------------------------------------------------------------//
//              end of file ResponseOperator.hh
//---------------------------------------------------------------------------//
