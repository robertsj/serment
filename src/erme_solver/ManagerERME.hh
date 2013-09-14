//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ManagerERME.hh
 *  @brief ManagerERME class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef erme_utils_MANAGERERME_HH_
#define erme_utils_MANAGERERME_HH_

//----------------------------------------------------------------------------//
/** @mainpage Serment: An Eigenvalue Response Matrix Code
 *
 *  @section sec_introduction Introduction
 *
 *  The eigenvalue response matrix method (ERMM) is a numerical approach
 *  for static reactor analysis.  The basic idea behind the method is
 *  to decompose a global domain in space and link the resulting independent
 *  \e nodes via approximate boundary conditions.  Because the nodes are
 *  completely independent in a computational sense, the "physical"
 *  decomposition results in natural computational domain decomposition
 *  for parallel computation.
 *
 *  Since nodes communicate only at the boundaries,
 *  a variety of transport methods, both deterministic and stochastic, can
 *  be used to generate the boundary conditions as long as a consistent
 *  approximation is used (<em>e.g.</em> a P1 approximation, etc.).
 *  Furthermore, different approximations can be used for different nodes,
 *  and varying levels of boundary approximation can be used throughout
 *  the domain.
 *
 *  A more detailed discussion of the methods used in Serment can be found
 *  in the theory documentation.  Details on implementation specifics can
 *  be found throughout the rest of this documentation.
 *
 *  @section install_sec Installation
 *
 *  @subsection step1 Step 1: Opening the box
 *
 * etc...
 */
//----------------------------------------------------------------------------//

#include "erme/StateERME.hh"
#include "erme/ResponseContainer.hh"
#include "erme_geometry/NodeList.hh"
#include "erme_geometry/NodePartitioner.hh"
#include "erme_response/ResponseServer.hh"
#include "erme_response/ResponseIndexer.hh"
#include "erme_solver/GlobalSolverBase.hh"
#include "utilities/DBC.hh"
#include "utilities/InputDB.hh"
#include "utilities/SP.hh"

/**
 *  @brief Namespace for higher level routines for organizing response
 *         matrix problems and their solutions.
 */
namespace erme_solver
{

/**
 *  @class ManagerERME
 *  @brief Manages solution of eigenvalue response matrix problem.
 */
class ManagerERME
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<ManagerERME>           	SP_manager;
  typedef detran_utilities::InputDB::SP_input         	SP_db;
  typedef erme::StateERME::SP_state                   	SP_state;
  typedef erme::ResponseContainer::SP_responsecontainer SP_responsecontainer;
  typedef erme_geometry::NodeList::SP_nodelist        	SP_nodelist;
  typedef erme_geometry::NodePartitioner              	Partitioner;
  typedef erme_response::ResponseServer::SP_server    	SP_server;
  typedef erme_response::ResponseIndexer::SP_indexer  	SP_indexer;
  typedef erme_solver::GlobalSolverBase::SP_solver    	SP_solver;
  typedef detran_utilities::vec_dbl                   	vec_dbl;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /// Constructor
  ManagerERME(int argc, char *argv[]);


  /// SP constructor
  static SP_manager Create(int argc, char *argv[]);

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /**
   *  @brief Construct the communicators
   *  @param db   User parameter db
   */
  void build_comm(SP_db db);

  /**
   *  @brief Construct the problem
   *  @param nodes  List of nodes
   */
  void build_erme(SP_nodelist nodes);

  /**
   *  @brief Solve the problem
   *
   *  The solver is built on-the-fly, which allows one to quickly
   *  change parameters via the input database for successive
   *  solves.
   */
  void solve();

  /// Getters
  //@{
  SP_indexer indexer() const { return d_indexer; }
  double get_keff() const { return d_state->k(); }
  double get_lambda() const { return d_state->lambda(); }
  vec_dbl get_residual_norms() const { return d_solver->residual_norms(); }
  SP_responsecontainer get_responses() const { return d_responses; }
  SP_server get_server() const { return d_server; }
  SP_state get_state() {return d_state;}
  SP_nodelist get_nodes() {return d_nodes;}
  SP_indexer get_indexer() {return d_indexer;}
  //@}

  /// Close libraries, etc.
  void finalize();

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Number of command line arguments
  int d_argc;
  /// Command line arguments
  char** d_argv;
  /// Parameter database
  SP_db d_db;
  /// Node list
  SP_nodelist d_nodes;
  /// Indexer
  SP_indexer d_indexer;
  /// Response server
  SP_server d_server;
  /// State vector
  SP_state d_state;
  /// Response container
  SP_responsecontainer d_responses;
  /// Global solver
  SP_solver d_solver;
  /// Is it built?
  bool d_is_built;

};

} // end namespace erme_solver

#endif // erme_utils_MANAGERERME_HH_

//----------------------------------------------------------------------------//
//              end of file ManagerERME.hh
//----------------------------------------------------------------------------//
