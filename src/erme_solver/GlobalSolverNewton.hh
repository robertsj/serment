//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  GlobalSolverNewton.hh
 *  @brief GlobalSolverNewton class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef erme_solver_GLOBALSOLVERNEWTON_HH_
#define erme_solver_GLOBALSOLVERNEWTON_HH_

#include "GlobalSolverBase.hh"
#include "Jacobian.hh"
#include "linear_algebra/NonlinearSolver.hh"

namespace erme_solver
{

/**
 *  @class GlobalSolverNewton
 *  @brief Solves the problem using Newton's method
 *
 *  Nonlinear problems are discussed briefly in @ref NonlinearSolver.  See
 *  @ref NonlinearResidual and @ref Jacobian for specifics on solving
 *  the eigenvalue response matrix equations via Newton's method.
 */
class GlobalSolverNewton: public GlobalSolverBase
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<GlobalSolverNewton>    SP_solver;
  typedef erme_solver::GlobalSolverBase               Base;
  typedef Base::SP_solver                             SP_base;
  typedef Jacobian::SP_jacobian                       SP_jacobian;
  typedef linear_algebra::NonlinearSolver::SP_solver  SP_nonlinear_solver;

  //--------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param db         Pointer to parameter database
   *  @param indexer    Pointer to response indexer
   *  @param server     Pointer to response server
   *  @param state      Pointer to state vector
   *  @param responses  Container of the responses
   */
  GlobalSolverNewton(SP_db                db,
                     SP_indexer           indexer,
                     SP_server            server,
                     SP_state             state,
                     SP_responsecontainer responses);

  /// Virtual destructor
  virtual ~GlobalSolverNewton(){}

  /// Solve
  void solve();

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Jacobian
  SP_jacobian d_jacobian;
  /// Jacobian preconditioner
  SP_jacobian d_preconditioner;
  /// Nonlinear solver
  SP_nonlinear_solver d_solver;
  /// Solver for initial guess
  SP_base d_initial;

};


} // end namespace erme_solver

#endif // erme_solver_GLOBALSOLVERNEWTON_HH_

//----------------------------------------------------------------------------//
//              end of file GlobalSolverNewton.hh
//----------------------------------------------------------------------------//
