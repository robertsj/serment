//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ResponseServer.hh
 *  @brief ResponseServer class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef erme_response_RESPONSESERVER_HH_
#define erme_response_RESPONSESERVER_HH_

#include "ResponseIndexer.hh"
#include "NodeResponse.hh"
#include "ResponseSource.hh"
#include "ResponseDatabase.hh"
#include "erme_geometry/NodeList.hh"
#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include "utilities/SP.hh"
#include <vector>
#include <string>

namespace erme_response
{

/**
 *  @class ResponseServer
 *  @brief Serve nodal responses to clients
 *
 *  A ResponseServer lives on a local communicator.  A server is in charge
 *  of one or more Node objects.  The nodal responses are produced by a
 *  ResponseSource that solves the local problems.  There is one
 *  source for each unique node.
 *
 */
/**
 *  @example erme_response/test/test_ResponseServer.cc
 *
 *  Test of ResponseServer class
 */
class ResponseServer
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<ResponseServer>  SP_server;
  typedef NodeResponse::SP_response             SP_response;
  typedef std::vector<SP_response>              vec_response;
  typedef ResponseSource::SP_source             SP_source;
  typedef std::vector<SP_source>                vec_source;
  typedef ResponseIndexer::SP_indexer           SP_indexer;
  typedef erme_geometry::NodeList::SP_nodelist  SP_nodelist;
  typedef ResponseDatabase::SP_rfdb             SP_rfdb;
  typedef detran_utilities::size_t              size_t;

  //--------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param nodes    Pointer to node list
   *  @param indexer  Pointer to indexer
   *  @param dbname   Filename of response database (optional)
   *  @param dborder  Interpolation order for database (optional)
   *  @param dbtype   Interpolation type for database (optional)
   */
  ResponseServer(SP_nodelist nodes,
                 SP_indexer  indexer,
                 std::string dbname = "",
                 size_t      dborder = 1,
                 std::string dbtype = "i");

  /// Update the eigenvalue and compute the new responses
  bool update(const double keff, int msg = 0);

  /**
   *  @brief Return a nodal response
   *  @param node   Local node index
   */
  SP_response response(size_t node);

  /// Ask server if it has updated responses (before reconstructing operators)
  bool is_updated() const;

  size_t number_unique_evaluations() const;
  size_t number_swapped_evaluations() const;
  /// Return the total time elapsed generating responses
  double response_time() const;

  /// Get the last keff not equal to the argument
  double last_keff(const double k) const;

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Nodes
  SP_nodelist d_nodes;
  /// Indexer
  SP_indexer d_indexer;
  /// Sources [size = number of unique local nodes]
  vec_source d_sources;
  /// Node response functions [size = number of unique local nodes]
  vec_response d_responses;
  /// Previously-computed response functions (if k is repeated)
  vec_response d_responses_1;
  /// Response database
  SP_rfdb d_rfdb;
  /// Current k-eigenvalue
  double d_keff;
  /// Last k-eigenvalue
  double d_keff_1;
  /// Flag indicating the server has new values
  bool d_is_updated;
  /// Number of unique evaluations (i.e. new keff values)
  size_t d_number_unique_evaluations;
  /// Number of evaluations that require updated operators but are in memory
  size_t d_number_swapped_evaluations;
  /// Response generation time
  double d_response_time;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  /**
   *  @brief Compute responses by dividing work explicitly among processes
   *
   *  Each response is treated the same, when in fact some responses (or
   *  nodes) might be inherently expensive.  A reduction is used to get
   *  the data, which probably is inefficient for lots of nodes and/or
   *  lots of response data.
   *
   */
  void update_explicit_work_share();

  /**
   *  @brief Compute responses by using a self-scheduling master-slave
   *         approach
   */
  void update_master_slave()
  {
    THROW("MASTER-SLAVE NOT IMPLEMENTED");
  }

};

} // end namespace erme_response

//----------------------------------------------------------------------------//
// INLINE MEMBER DEFINITIONS
//----------------------------------------------------------------------------//

#include "ResponseServer.i.hh"

#endif // erme_response_RESPONSESERVER_HH_

//----------------------------------------------------------------------------//
//              end of file ResponseServer.hh
//----------------------------------------------------------------------------//
