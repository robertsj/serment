//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  GlobalSolverSteffensen.hh
 *  @brief GlobalSolverSteffensen class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef erme_solver_GLOBALSOLVERSTEFFENSEN_HH_
#define erme_solver_GLOBALSOLVERSTEFFENSEN_HH_

#include "GlobalSolverBase.hh"

namespace erme_solver
{

/**
 *  @class GlobalSolverPicard
 *  @brief Solves the problem using Steffensen  iteration
 *
 *  Steffensen's method solves the fixed point problem @f$ x' = f(x) @f$
 *  via the modified process
 *  @f[
 *      x' = g(x) = x - (x-f(x))^2 / (f(f(x))-2*f(x)-x) \, .
 *  @f]
 *  This can be shown to be second order.
 *
 *  Note, an @ref SteffensenUpdate is an Aitken-based versio for
 *  accelerating @ref GlobalSolverPicard.  The implementation here should
 *  in theory given identical results.
 */
class GlobalSolverSteffensen: public GlobalSolverBase
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<GlobalSolverSteffensen>  SP_solver;
  typedef GlobalSolverBase                              Base;

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
  GlobalSolverSteffensen(SP_db                db,
                         SP_indexer           indexer,
                         SP_server            server,
                         SP_state             state,
                         SP_responsecontainer responses);

  /// Virtual destructor
  virtual ~GlobalSolverSteffensen(){}

  /// Solve
  void solve();

private:

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  void display(const size_t it,
               const double norm,
               const double keff,
               const double lambda);

};

} // end namespace erme_solver

#endif /* erme_solver_GLOBALSOLVERSTEFFENSEN_HH_ */

//----------------------------------------------------------------------------//
//              end of file GlobalSolverSteffensen.hh
//----------------------------------------------------------------------------//
