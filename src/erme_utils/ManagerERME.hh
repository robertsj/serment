//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ManagerERME.hh
 *  @brief  ManagerERME class definition
 *  @author Jeremy Roberts
 *  @date   Sep 1, 2012
 */
//---------------------------------------------------------------------------//

#ifndef erme_utils_MANAGERERME_HH_
#define erme_utils_MANAGERERME_HH_

//---------------------------------------------------------------------------//
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
//---------------------------------------------------------------------------//

#include "erme/StateERME.hh"
#include "erme/ResponseMatrix.hh"
#include "erme/Connect.hh"
#include "erme/FissionOperator.hh"
#include "erme/AbsorptionOperator.hh"
#include "erme/LeakageOperator.hh"
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
 *         matrix problems and their solutions
 */
namespace erme_utils
{

/**
 *  @class ManagerERME
 *  @brief Manages solution of eigenvalue response matrix problem.
 */
class ManagerERME
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<ManagerERME>           SP_manager;
  typedef detran_utilities::InputDB::SP_input         SP_db;
  typedef erme::StateERME::SP_state                   SP_state;
  typedef erme::ResponseMatrix::SP_responsematrix     SP_R;
  typedef erme::Connect::SP_connect                   SP_M;
  typedef erme::FissionOperator::SP_fission           SP_F;
  typedef erme::AbsorptionOperator::SP_absorption     SP_A;
  typedef erme::LeakageOperator::SP_leakage           SP_L;
  typedef erme_geometry::NodeList::SP_nodelist        SP_nodelist;
  typedef erme_geometry::NodePartitioner              Partitioner;
  typedef erme_response::ResponseServer::SP_server    SP_server;
  typedef erme_response::ResponseIndexer::SP_indexer  SP_indexer;
  typedef erme_solver::GlobalSolverBase::SP_solver    SP_solver;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /// Constructor
  ManagerERME(int argc, char *argv[], SP_db db);

  /// SP constructor
  static SP_manager Create(int argc, char *argv[], SP_db db);

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

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

  /// Return the indexer
  SP_indexer indexer() const { return d_indexer; }

  /// Close libraries, etc.
  void finalize();

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

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
  /// Response matrix
  SP_R d_R;
  /// Connectivity matrix
  SP_M d_M;
  /// Fission operator
  SP_F d_F;
  /// Absorption operator
  SP_A d_A;
  /// Leakage operator
  SP_L d_L;
  /// Global solver
  SP_solver d_solver;


};

} // end namespace erme_utils

#endif // erme_utils_MANAGERERME_HH_

//---------------------------------------------------------------------------//
//              end of file ManagerERME.hh
//---------------------------------------------------------------------------//
