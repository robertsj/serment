//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   GlobalSolverBase.hh
 * \brief  GlobalSolverBase class definition
 * \author Jeremy Roberts
 * \date   Sep 4, 2012
 */
//---------------------------------------------------------------------------//

#ifndef GLOBALSOLVERBASE_HH_
#define GLOBALSOLVERBASE_HH_

#include "erme_response/ResponseServer.hh"
#include "erme/StateERME.hh"
#include "erme/ResponseMatrix.hh"
#include "erme/Connect.hh"
#include "erme/FissionOperator.hh"
#include "erme/AbsorptionOperator.hh"
#include "erme/LeakageOperator.hh"
#include "utilities/DBC.hh"
#include "utilities/SP.hh"

/**
 *  @namespace erme_solver
 *  @brief Contains the solver package
 */
namespace erme_solver
{

/*!
 *  \class GlobalSolverBase
 */
class GlobalSolverBase
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::InputDB::SP_input       SP_db;
  typedef detran_utilities::SP<GlobalSolverBase>    SP_solver;
  typedef erme_response::ResponseServer::SP_server  SP_server;
  typedef erme::StateERME::SP_state                 SP_state;
  typedef erme::ResponseMatrix::SP_responsematrix   SP_R;
  typedef erme::Connect::SP_connect                 SP_M;
  typedef erme::FissionOperator::SP_fission         SP_F;
  typedef erme::AbsorptionOperator::SP_absorption   SP_A;
  typedef erme::LeakageOperator::SP_leakage         SP_L;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   *  \param db     Pointer to parameter database
   *  \param server Pointer to response server
   *  \param state  Pointer to state vector
   *  \param R      Pointer to response matrix
   *  \param M      Pointer to connectivity matrix
   *  \param F      Pointer to fission operator
   *  \param A      Pointer to absorption operator
   *  \param L      Pointer to leakage operator
   */
  GlobalSolverBase(SP_db db, SP_server server, SP_state state,
                   SP_R R, SP_M M, SP_F F, SP_A A, SP_L L);

  /// Solve
  virtual void solve() = 0;

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

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

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//


};


} // end namespace erme_solver

#endif // GLOBALSOLVERBASE_HH_ 

//---------------------------------------------------------------------------//
//              end of file GlobalSolverBase.hh
//---------------------------------------------------------------------------//
