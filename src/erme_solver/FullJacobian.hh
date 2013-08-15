//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  FullJacobian.hh
 *  @brief FullJacobian class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef erme_solver_FULLJACOBIAN_HH_
#define erme_solver_FULLJACOBIAN_HH_

#include "OperatorMR.hh"
#include "linear_algebra/JacobianBase.hh"
#include "linear_algebra/MatrixShell.hh"
#include "erme_response/ResponseIndexer.hh"
#include "erme_response/ResponseServer.hh"
#include "erme/StateERME.hh"
#include "erme/ResponseMatrix.hh"
#include "erme/Connect.hh"
#include "erme/FissionOperator.hh"
#include "erme/AbsorptionOperator.hh"
#include "erme/LeakageOperator.hh"
#include "erme/ResponseContainer.hh"
#include "utilities/DBC.hh"
#include "utilities/SP.hh"
#include "utilities/Definitions.hh"

namespace erme_solver
{

//----------------------------------------------------------------------------//
/**
 *  @class FullJacobian
 *  @brief This constructs a full, possibly approximate Jacobian
 *
 *  This class should primarily be used to provide an accurate preconditioner
 *  for use with the shell Jacobian.  The finite difference component is
 *  optional left out to avoid extra response evaluations
 *
 *  For a more detailed description of the Jacobian operator, see
 *  @ref Jacobian
 */
class FullJacobian: public linear_algebra::JacobianBase
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef erme_response::ResponseServer::SP_server        SP_server;
  typedef erme::ResponseContainer::SP_responsecontainer   SP_responsecontainer;
  typedef erme::ResponseMatrix::SP_responsematrix         SP_R;
  typedef erme::Connect::SP_connect                       SP_M;
  typedef erme::FissionOperator::SP_fission               SP_F;
  typedef erme::AbsorptionOperator::SP_absorption         SP_A;
  typedef erme::LeakageOperator::SP_leakage               SP_L;
  typedef linear_algebra::Vector                          Vector;
  typedef Vector::SP_vector                               SP_vector;
  typedef detran_utilities::vec_dbl                       vec_dbl;
  typedef detran_utilities::vec_int                       vec_int;
  typedef detran_utilities::size_t                        size_t;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  FullJacobian(SP_server              server,
               SP_responsecontainer   responses,
               const double           eps = 1.0e-8);

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// Update the Jacobian
  void update(SP_vector x);

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Response server
  SP_server d_server;
  //@{
  /// Operators
  SP_R d_R;
  SP_M d_M;
  SP_F d_F;
  SP_A d_A;
  SP_L d_L;
  //@}
  /// Finite difference epsilon
  double d_eps;
  /// Unknown vector
  SP_vector d_x;
  /// Eigenvalues
  double d_k;
  double d_lambda;
  /// Finite difference of (d/dk)[(MR) * J]
  SP_vector d_fd_MR;
  /// Finite difference of (d/dk)[(F - k*L) * J]
  double d_fd_FAL;
  /// Full and truncated local sizes
  size_t d_m_full;
  size_t d_m;

};

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of file FullJacobian.hh
//----------------------------------------------------------------------------//




#endif /* FULLJACOBIAN_HH_ */