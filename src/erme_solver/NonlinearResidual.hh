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
#include "erme/ResponseContainer.hh"
#include "linear_algebra/Vector.hh"
#include "linear_algebra/NonlinearResidualBase.hh"
#include "utilities/SP.hh"

namespace erme_solver
{

/**
 *  @class NonlinearResidual
 *  @brief Provides the residual of the response matrix equations
 *
 *  The nonlinear residual is defined as
 *  @f[
 *    \mathbf{f(x)} =
 *      \left [\begin{array}{c}
 *              (\mathbf{M}\mathbf{R}(k)-\lambda \mathbf{I}) \mathbf{J_-} \\
 *              \mathbf{F}(k)\mathbf{J_-} - (k\mathbf{L}(k)\mathbf{J_-} ) \\
 *              \frac{1}{2} \mathbf{J^T_-} \mathbf{J_-} - \frac{1}{2}
 *              \end{array}
 *       \right ]  = \mathbf{0} \, ,
 *  @f]
 *  which is used to assess convergence of all solvers and as part of
 *  Newton schemes.
 */
class NonlinearResidual: public linear_algebra::NonlinearResidualBase
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<NonlinearResidual>         SP_residual;
  typedef erme_response::ResponseServer::SP_server        SP_server;
  typedef erme::ResponseContainer::SP_responsecontainer   SP_responses;
  typedef linear_algebra::Vector                          Vector;
  typedef linear_algebra::Vector::SP_vector               SP_vector;
  typedef OperatorMR::SP_matrix                           SP_matrix;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param    server      Response server
   *  @param    responses   Container of responses
   */
  NonlinearResidual(SP_server server, SP_responses responses);

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// Evaluate the residual vector given a vector of unknown
  void evaluate(Vector *x, Vector *f);

  /// Computes \f$ L_2 \f$ of the residual given a vector of all unknowns
  double compute_norm(Vector *x);

  /// Computes \f$ L_2 \f$ of the residual given current and eigenvalues
  double compute_norm(Vector *J, const double k, const double l);

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Response server
  SP_server d_server;
  /// Response container
  SP_responses d_responses;
  /// Response matrix operator
  SP_matrix d_MR;

  /// Update the response operators given a new eigenvalue
  void update_response(const double keff);
};

} // end namespace erme_solver

#endif // erme_solver_NONLINEARRESIDUAL_HH_

//----------------------------------------------------------------------------//
//              end of file NonlinearResidual.hh
//----------------------------------------------------------------------------//
