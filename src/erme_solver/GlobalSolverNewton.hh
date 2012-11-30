//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   GlobalSolverNewton.hh
 *  @brief  GlobalSolverNewton
 *  @author Jeremy Roberts
 *  @date   Nov 29, 2012
 */
//---------------------------------------------------------------------------//

#ifndef erme_solver_GLOBALSOLVERNEWTON_HH_
#define erme_solver_GLOBALSOLVERNEWTON_HH_

#include "GlobalSolverBase.hh"
#include "OperatorMR.hh"

namespace erme_solver
{

/**
 *  @class GlobalSolverNewton
 *  @brief Solves the problem using Picard (fixed point) iteration
 *
 */
class GlobalSolverNewton: public GlobalSolverBase
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<GlobalSolverNewton>    SP_solver;
  typedef GlobalSolverBase                            Base;
  typedef erme

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param server Pointer to response server
   *  @param state  Pointer to state vector
   *  @param R      Pointer to response matrix
   *  @param M      Pointer to connectivity matrix
   *  @param F      Pointer to fission operator
   *  @param A      Pointer to absorption operator
   *  @param L      Pointer to leakage operator
   */
  GlobalSolverNewton(SP_db db, SP_indexer indexer, SP_server server,
                     SP_state state,
                     SP_R R, SP_M M, SP_F F, SP_A A, SP_L L);

  /// Virtual destructor
  virtual ~GlobalSolverNewton(){}

  /// Solve
  void solve();

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Jacobian
  SP_Jacobian d_jacobian;
  /// MR operator
  SP_MR d_MR;
  /// Temporary working vectors
  SP_vector d_J0;
  SP_vector d_J1;

};


} // end namespace erme_solver

#endif // erme_solver_GLOBALSOLVERNEWTON_HH_

//---------------------------------------------------------------------------//
//              end of file GlobalSolverNewton.hh
//---------------------------------------------------------------------------//
