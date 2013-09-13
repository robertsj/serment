//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Jacobian.hh
 *  @brief Jacobian class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef erme_solver_JACOBIAN_HH_
#define erme_solver_JACOBIAN_HH_

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
 *  @class Jacobian
 *  @brief This provides the action of the Jacobian
 *         (\f$ \mathbf{f}' \f$ )on a vector \f$ \mathbf{f} \f$
 *
 *  The resulting vector is \f$ \mathbf{f}'\mathbf{f} \f$, or "fpf" where "p"
 *  denotes the prime.
 *
 *  Note, the vector \f$ \vec{f} \f$  does not necessarily begin as the
 *  residual but rather as a scaled residual. The PETSc manual does not seem
 *  to state this explicitly, but a little trial and error confirmed this---so
 *  don't fret when debugging!
 *
 *  The Jacobian is defined as
 *  \f[
 *      \mathbf{f'(x)} =
 *        \left [\begin{array}{ccc}
 *          (\mathbf{M}\mathbf{R}-\lambda \mathbf{I})              &
 *              \mathbf{M}\mathbf{R_k}\mathbf{J_-}                     &
 *                  \mathbf{J_-}                                           \\
 *          (\mathbf{F}-k\mathbf{L})                               &
 *              (\mathbf{F_k}-k\mathbf{L_k}-\mathbf{L}) \mathbf{J_-}   &
 *                  0                                                      \\
 *          \mathbf{J^T_-}                                         &
 *              0                                                      &
 *                 0
 *        \end{array} \right ]  \, .
 *      \label{eq:jacobian}
 *  \f]
 *  For \f$ \mathbf{R}(k) \f$ of size \f$ m\times m \f$, the Jacobian is
 *  of size \f$ (m+2)\times(m+2) \f$.
 *
 *  Most of the Jacobian is known <EM> a priori</EM>, and so most of the
 *  action can be computed directly.  Only the first \f$ m-1 \f$ rows of the
 *  \f$ (m-1) \f$th column require approximations via finite differences.
 *
 */
class Jacobian: public linear_algebra::JacobianBase
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef erme_response::ResponseServer::SP_server        SP_server;
  typedef erme_response::ResponseIndexer::SP_indexer      SP_indexer;
  typedef erme::ResponseContainer::SP_responsecontainer   SP_responsecontainer;
  typedef erme::ResponseMatrix::SP_responsematrix         SP_R;
  typedef erme::Connect::SP_connect                       SP_M;
  typedef erme::FissionOperator::SP_fission               SP_F;
  typedef erme::AbsorptionOperator::SP_absorption         SP_A;
  typedef erme::LeakageOperator::SP_leakage               SP_L;
  typedef OperatorMR::SP_MR                               SP_MR;
  typedef linear_algebra::Vector                          Vector;
  typedef Vector::SP_vector                               SP_vector;
  typedef detran_utilities::vec_dbl                       vec_dbl;
  typedef detran_utilities::vec_int                       vec_int;
  typedef detran_utilities::size_t                        size_t;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param    server      response server
   *  @param    indexer     response indexer
   *  @param    responses   container of response operators
   *  @param    eps         delta-k for finite difference
   */
  Jacobian(SP_server              server,
           SP_indexer             indexer,
           SP_responsecontainer   responses,
           const double           eps = 1.0e-8);

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  void multiply(Vector &v_in, Vector &v_out);
  void multiply_transpose(Vector &v_in, Vector &v_out);

  /// Update the Jacobian
  void update(SP_vector x);

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Response server
  SP_server d_server;
  /// Response indexer
  SP_indexer d_indexer;
  //@{
  /// Operators
  SP_R d_R;
  SP_M d_M;
  SP_F d_F;
  SP_A d_A;
  SP_L d_L;
  SP_MR d_MR;
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

  /// Shell matrix for Jacobian action
  class Shell: public linear_algebra::MatrixShell
  {
  public:
    Shell(const size_type m, Jacobian &J)
      : linear_algebra::MatrixShell(m, m, this)
      , d_J(J)
    {
      /* ... */
    }
    void multiply(Vector &v_in, Vector &v_out)
    {
      d_J.multiply(v_in, v_out);
    }
    void multiply_transpose(Vector &v_in, Vector &v_out)
    {
      d_J.multiply_transpose(v_in, v_out);
    }
  private:
    Jacobian& d_J;
  };

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  /// Update the response operators given a new eigenvalue
  void update_response(const double keff);

};

} // end namespace detran

#endif // erme_solver_JACOBIAN_HH_

//----------------------------------------------------------------------------//
//              end of file Jacobian.hh
//----------------------------------------------------------------------------//
