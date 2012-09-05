//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   GlobalSolverPicard.hh
 * \brief  GlobalSolverPicard class definition
 * \author Jeremy Roberts
 * \date   Sep 4, 2012
 */
//---------------------------------------------------------------------------//

#ifndef GLOBALSOLVERPICARD_HH_
#define GLOBALSOLVERPICARD_HH_

#include "GlobalSolverBase.hh"

namespace erme_solver
{

class GlobalSolverPicard: public GlobalSolverBase
{

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran::SP<GlobalSolverPicard>              SP_solver;
  typedef GlobalSolverBase                            Base;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   *  \param server Pointer to response server
   *  \param state  Pointer to state vector
   *  \param R      Pointer to response matrix
   *  \param M      Pointer to connectivity matrix
   *  \param F      Pointer to fission operator
   *  \param A      Pointer to absorption operator
   *  \param L      Pointer to leakage operator
   */
  GlobalSolverPicard(SP_db db, SP_server server, SP_state state,
                     SP_R R, SP_M M, SP_F F, SP_A A, SP_L L);

  /// Solve
  virtual void solve();



};


} // end namespace erme_solver

#endif // GLOBALSOLVERPICARD_HH_ 

//---------------------------------------------------------------------------//
//              end of file GlobalSolverPicard.hh
//---------------------------------------------------------------------------//
