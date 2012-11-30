//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   GlobalSolverBase.hh
 *  @brief  GlobalSolverBase class definition
 *  @author Jeremy Roberts
 *  @date   Sep 4, 2012
 */
//---------------------------------------------------------------------------//

#ifndef erme_solver_GLOBALSOLVERBASE_HH_
#define erme_solver_GLOBALSOLVERBASE_HH_

#include "NonlinearResidual.hh"
#include "erme_response/ResponseIndexer.hh"
#include "erme_response/ResponseServer.hh"
#include "erme/StateERME.hh"
#include "erme/ResponseMatrix.hh"
#include "erme/Connect.hh"
#include "erme/FissionOperator.hh"
#include "erme/AbsorptionOperator.hh"
#include "erme/LeakageOperator.hh"
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

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<GlobalSolverBase>      SP_solver;
  typedef detran_utilities::InputDB::SP_input         SP_db;
  typedef erme_response::ResponseIndexer::SP_indexer  SP_indexer;
  typedef erme_response::ResponseServer::SP_server    SP_server;
  typedef erme::StateERME::SP_state                   SP_state;
  typedef erme::ResponseMatrix::SP_responsematrix     SP_R;
  typedef erme::Connect::SP_connect                   SP_M;
  typedef erme::FissionOperator::SP_fission           SP_F;
  typedef erme::AbsorptionOperator::SP_absorption     SP_A;
  typedef erme::LeakageOperator::SP_leakage           SP_L;
  typedef OperatorMR::SP_MR                           SP_MR;
  typedef NonlinearResidual::SP_residual              SP_residual;
  typedef linear_algebra::Vector::SP_vector           SP_vector;
  typedef detran_utilities::vec_dbl                   vec_dbl;
  typedef detran_utilities::vec_int                   vec_int;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param db       Pointer to parameter database
   *  @param indexer  Pointer to response indexer
   *  @param server   Pointer to response server
   *  @param state    Pointer to state vector
   *  @param R        Pointer to response matrix
   *  @param M        Pointer to connectivity matrix
   *  @param F        Pointer to fission operator
   *  @param A        Pointer to absorption operator
   *  @param L        Pointer to leakage operator
   */
  GlobalSolverBase(SP_db db, SP_indexer indexer, SP_server server,
                   SP_state state,
                   SP_R R, SP_M M, SP_F F, SP_A A, SP_L L);

  /// Pure virtual destructor
  virtual ~GlobalSolverBase() = 0;

  /// Solve
  virtual void solve() = 0;

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Parameter database
  SP_db d_db;
  /// Response index
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
  /// Residual norm computer
  SP_residual d_residual;
  /// Maximum (outer) iterations
  int d_maximum_iterations;
  /// Convergence tolerance
  double d_tolerance;
  /// Residual norm for each (outer) iteration
  vec_dbl d_residual_norms;
  /// Number of inners for each (outer) iteration
  vec_int d_inner_iterations;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /// Update the response operators given a new eigenvalue
  void update_response(const double keff)
  {
    /// Give the server the new keff
    d_server->update(keff);
    /// Fill the operators with updated values
    d_R->update();
    d_F->update();
    d_A->update();
    d_L->update();
  }

  /// Monitor the convergence
  void monitor();

};

} // end namespace erme_solver

#endif // erme_solver_GLOBALSOLVERBASE_HH_

//---------------------------------------------------------------------------//
//              end of file GlobalSolverBase.hh
//---------------------------------------------------------------------------//
