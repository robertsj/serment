//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Jacobian.hh
 *  @brief  Jacobian
 *  @author Jeremy Roberts
 *  @date   Nov 29, 2012
 */
//---------------------------------------------------------------------------//

#ifndef erme_solver_JACOBIAN_HH_
#define erme_solver_JACOBIAN_HH_

#include "linear_algebra/MatrixShell.hh"
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

namespace erme_solver
{

//---------------------------------------------------------------------------//
/**
 *  @class Jacobian
 *  @brief This performs the action of the Jacobian
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
class Jacobian: public linear_algebra::MatrixShell
{

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef erme::StateERME::SP_state                   SP_state;
  typedef erme::ResponseMatrix::SP_responsematrix     SP_R;
  typedef erme::Connect::SP_connect                   SP_M;
  typedef erme::FissionOperator::SP_fission           SP_F;
  typedef erme::AbsorptionOperator::SP_absorption     SP_A;
  typedef erme::LeakageOperator::SP_leakage           SP_L;
  typedef OperatorMR::SP_MR                           SP_MR;
  typedef linear_algebra::Vector::SP_vector           SP_vector;
  typedef detran_utilities::vec_dbl                   vec_dbl;
  typedef detran_utilities::vec_int                   vec_int;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  Jacobian(SP_state state, SP_R R, SP_M M, SP_F F, SP_A A, SP_L L);

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /**
   *  @brief Matrix-vector multiplication
   *  @param x  Input vector
   *  @param y  Output vector
   */
  PetscErrorCode shell_multiply(Vec x, Vec y);

  /**
   *  @brief Matrix-vector multiplication using matrix transpose.
   *  @param x  Input vector
   *  @param y  Output vector
   */
  PetscErrorCode shell_multiply_transpose(Vec x, Vec y)
  {
    THROW("TRANSPOSE JACOBIAN UNAVAILABLE");
    return 0;
  }

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// MR operator
  SP_MR d_MR;
  /// Temporary working vectors
  SP_vector d_J0;
  SP_vector d_J1;

};


} // end namespace detran

#endif // erme_solver_JACOBIAN_HH_

//---------------------------------------------------------------------------//
//              end of file Jacobian.hh
//---------------------------------------------------------------------------//
