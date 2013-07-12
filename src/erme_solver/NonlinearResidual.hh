//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  NonlinearResidual.hh
 *  @brief NonlinearResidual
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef erme_solver_NONLINEARRESIDUAL_HH_
#define erme_solver_NONLINEARRESIDUAL_HH_

#include "OperatorMR.hh"
#include "erme/ResponseMatrix.hh"
#include "erme/Connect.hh"
#include "erme/FissionOperator.hh"
#include "erme/AbsorptionOperator.hh"
#include "erme/LeakageOperator.hh"
#include "linear_algebra/Vector.hh"
#include "linear_algebra/NonlinearResidualBase.hh"
#include "utilities/SP.hh"

namespace erme_solver
{

class NonlinearResidual: public linear_algebra::NonlinearResidualBase
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<NonlinearResidual>     SP_residual;
  typedef erme::ResponseMatrix::SP_responsematrix     SP_R;
  typedef erme::Connect::SP_connect                   SP_M;
  typedef erme::FissionOperator::SP_fission           SP_F;
  typedef erme::AbsorptionOperator::SP_absorption     SP_A;
  typedef erme::LeakageOperator::SP_leakage           SP_L;
  typedef linear_algebra::Vector                      Vector;
  typedef linear_algebra::Vector::SP_vector           SP_vector;
  typedef OperatorMR::SP_matrix                       SP_matrix;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param R      Pointer to response matrix
   *  @param M      Pointer to connectivity matrix
   *  @param F      Pointer to fission operator
   *  @param A      Pointer to absorption operator
   *  @param L      Pointer to leakage operator
   */
  NonlinearResidual(SP_R R, SP_M M, SP_F F, SP_A A, SP_L L)
    : d_R(R), d_M(M), d_F(F), d_A(A), d_L(L)
  {
    d_MR = new OperatorMR(d_R, d_M);
  }

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /**
   *  @brief Computes the \f$ L_2 \f$ norm of the nonlinear residual.
   *
   *   The nonlinear residual is defined as
   *   @f[
   *     \mathbf{f(x)} =
   *       \left [\begin{array}{c}
   *              (\mathbf{M}\mathbf{R}(k)-\lambda \mathbf{I}) \mathbf{J_-} \\
   *              \mathbf{F}(k)\mathbf{J_-} - (k\mathbf{L}(k)\mathbf{J_-} ) \\
   *              \frac{1}{2} \mathbf{J^T_-} \mathbf{J_-} - \frac{1}{2}
   *              \end{array}
   *       \right ]  = \mathbf{0} \, ,
   *   @f]
   *   which is the same as used in the Newton-based schemes.  The \f$ L_2 \f$
   *   norm is then \f$ \sqrt{ \mathbf{f(x)}^T \mathbf{f(x)} } \f$.
   *
   *   @param x   vector of boundary unknowns with
   *              \f$ k \f$ and \f$ \lambda \f$
   */
  double compute_norm(Vector &x);

  /**
   *  @brief Computes the \f$ L_2 \f$ norm of the nonlinear residual.
   *
   *  This version offers an interface for Picard iteration.
   */
  double compute_norm(Vector &x, const double k, const double l);

  /// Evaluate the residual vector
  void evaluate(Vector &x, Vector &f);

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  SP_R d_R;
  SP_M d_M;
  SP_F d_F;
  SP_A d_A;
  SP_L d_L;
  SP_matrix d_MR;


};

} // end namespace erme_solver

#include "NonlinearResidual.i.hh"

#endif // erme_solver_NONLINEARRESIDUAL_HH_

//----------------------------------------------------------------------------//
//              end of file NonlinearResidual.hh
//----------------------------------------------------------------------------//
