//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ResponseSource.hh
 *  @brief ResponseSource class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef erme_response_RESPONSESOURCE_HH_
#define erme_response_RESPONSESOURCE_HH_

#include "NodeResponse.hh"
#include "ResponseIndex.hh"
#include "ResponseIndexer.hh"
#include "erme_geometry/Node.hh"
#include "utilities/DBC.hh"
#include "utilities/SP.hh"

namespace erme_response
{

/**
 *  @class ResponseSource
 *  @brief Abstract response source
 *
 *  A ResponseSource provides its ResponseServer with responses for use
 *  in the global solve.  Each ResponseSource is unique for a given Node.
 *  The ResponseSource represents the interface between Serment and local
 *  solvers such as Detran.  Each concrete Node implementation must have
 *  a corresponding concrete ResponseSource and
 *  ResponseSourceFactory::build specialization.
 *
 *
 */
class ResponseSource
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<ResponseSource>  SP_source;
  typedef unsigned int                          size_t;
  typedef NodeResponse::SP_response             SP_response;
  typedef ResponseIndexer::SP_indexer           SP_indexer;
  typedef erme_geometry::Node::SP_node          SP_node;

  //--------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param node   Pointer to node object for which this source generates
   *                responses
   */
  ResponseSource(SP_node node, SP_indexer indexer)
    : d_node(node)
    , d_indexer(indexer)
    , d_keff(1.0)
  {
    // Preconditions
    Require(node);
    Require(indexer);
  }

  /// Virtual destructor
  virtual ~ResponseSource(){}

  /// Update the k-eigenvalue
  void update(const double keff)
  {
    d_keff = keff;
  }

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE
  //--------------------------------------------------------------------------//

  /**
   *  @brief Compute a response for the requested incident index
   *
   *  The client passes the response to be updated and the index of
   *  the corresponding incident condition.  The client is then
   *  responsible for moving data from the sources to the server.
   *
   *  @param response   Pointer to response object to be updated
   *  @param index      Response indices
   */
  virtual void compute(SP_response response, const ResponseIndex &index) = 0;

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Node
  SP_node d_node;
  /// Response indexer
  SP_indexer d_indexer;
  /// K-eigenvalue
  double d_keff;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

};

} // end namespace erme_response

#endif // erme_response_RESPONSESOURCE_HH_

//----------------------------------------------------------------------------//
//              end of file ResponseSource.hh
//----------------------------------------------------------------------------//
