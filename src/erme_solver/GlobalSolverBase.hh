//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  GlobalSolverBase.hh
 *  @brief GlobalSolverBase class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef erme_solver_GLOBALSOLVERBASE_HH_
#define erme_solver_GLOBALSOLVERBASE_HH_

#include "OperatorMR.hh"
#include "NonlinearResidual.hh"
#include "erme_response/ResponseIndexer.hh"
#include "erme_response/ResponseServer.hh"
#include "erme/StateERME.hh"
#include "erme/ResponseContainer.hh"
#include "utilities/DBC.hh"
#include "utilities/SP.hh"
#include "utilities/Definitions.hh"

/**
 *  @namespace erme_solver
 *  @brief Contains solvers for eigenvalue response matrix problems
 */
namespace erme_solver
{

/**
 *  @class GlobalSolverBase
 *  @class Base class for eigenvalue response matrix solver
 */
class GlobalSolverBase
{

public:

  //--------------------------------------------------------------------------//
  // ENUMERATIONS
  //--------------------------------------------------------------------------//
//
//  enum SOLVER_STATUS
//  {
//    CONTINUE, COMPLETED, END_SOLVER_STATUS
//  };

  const static int CONTINUE = 1984;
  const static int COMPLETED = 2013;

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<GlobalSolverBase>      		SP_solver;
  typedef detran_utilities::InputDB::SP_input         		SP_db;
  typedef erme_response::ResponseIndexer::SP_indexer  		SP_indexer;
  typedef erme_response::ResponseServer::SP_server    		SP_server;
  typedef erme::StateERME::SP_state                   		SP_state;
  typedef erme::ResponseContainer::SP_responsecontainer		SP_responsecontainer;
  typedef erme::ResponseMatrix::SP_responsematrix     		SP_R;
  typedef erme::Connect::SP_connect                   		SP_M;
  typedef erme::FissionOperator::SP_fission           		SP_F;
  typedef erme::AbsorptionOperator::SP_absorption     		SP_A;
  typedef erme::LeakageOperator::SP_leakage           		SP_L;
  typedef OperatorMR::SP_MR                           		SP_MR;
  typedef NonlinearResidual::SP_residual              		SP_residual;
  typedef linear_algebra::Vector                          Vector;
  typedef Vector::SP_vector           		                SP_vector;
  typedef detran_utilities::vec_dbl                   		vec_dbl;
  typedef detran_utilities::vec_int                   		vec_int;
  typedef detran_utilities::vec_size_t                    vec_size_t;
  typedef detran_utilities::size_t                        size_t;

  //--------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param db       	Pointer to parameter database
   *  @param indexer  	Pointer to response indexer
   *  @param server   	Pointer to response server
   *  @param state    	Pointer to state vector
   *  @param responses	Container of the responses
   */
  GlobalSolverBase(SP_db 									db,
                   SP_indexer 						indexer,
                   SP_server 							server,
                   SP_state 							state,
                   SP_responsecontainer 	responses);

  /// Pure virtual destructor
  virtual ~GlobalSolverBase() = 0;

  /// Solve
  virtual void solve() = 0;

  /// Get the residuals
  vec_dbl residual_norms() const;


  //@{
  /// Get the number of inner or outer iterations
  size_t number_outer_iterations() const;
  size_t number_inner_iterations() const;
  vec_size_t number_inner_iterations_per_outer() const;
  //@}

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Parameter database
  SP_db d_db;
  /// Response index
  SP_indexer d_indexer;
  /// Response server
  SP_server d_server;
  /// State vector
  SP_state d_state;
  /// Container of all responses
  SP_responsecontainer d_responses;
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
  /// Residual norm computer
  SP_residual d_residual;
  /// Maximum (outer) iterations
  int d_maximum_iterations;
  /// Convergence tolerance
  double d_tolerance;
  /// Residual norm for each (outer) iteration
  vec_dbl d_residual_norms;
  /// Number of outer iterations
  size_t d_number_outer_iterations;
  /// Number of inner iterations
  size_t d_number_inner_iterations;
  /// Number of inners for each (outer) iteration
  vec_size_t d_number_inner_iterations_per_outer;
  /// Local size of unknown vector
  size_t d_local_size;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  /// Update the response operators given a new eigenvalue
  void update_response(const double keff);

  /// Monitor the convergence
  //void monitor();

  /// Set initial guess to be unity in zeroth order moments
  void setup_initial_current(Vector &x);

  /// Display responses
  void display_response(std::string s);

};

} // end namespace erme_solver

#endif // erme_solver_GLOBALSOLVERBASE_HH_

//----------------------------------------------------------------------------//
//              end of file GlobalSolverBase.hh
//----------------------------------------------------------------------------//
