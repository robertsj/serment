//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PowerIter.hh
 * \author Jeremy Roberts
 * \date   11/23/2010
 * \brief  PowerIter class definition.
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 177                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-12-09 16:31:49 -0500 (Fri, 09 Dec 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef POWERITER_HH
#define POWERITER_HH

#include "serment_config.h"

#include <iostream>
#include <cmath>
#include "petsc.h"

#include "LinAlg.hh"
#include "GlobalInput.hh"
#include "ResponseFunctionServer.hh"
#include "ResponseMatrix.hh"
#include "ResponseMatrixFull.hh"
#include "AbsorptionResponse.hh"
#include "FissionResponse.hh"
#include "LeakageResponse.hh"
#include "ConnectMatrix.hh"
#include "Connect2dCart.hh"
#include "GlobalProblem.hh"
#include "GlobalSolver.hh"

#include "InnerIterBase.hh"
#include "InnerIterPower.hh"

#ifdef SERMENT_ENABLE_SLEPC
#include "InnerIterSLEPc.hh"
#endif

using namespace std;

//===========================================================================//
/*!
 * \class PowerIter
 * \brief This class solves a \ref GlobalProblem via power iteration.
 *
 * The power (iteration) method is a standard procedure for finding the
 * largest eigenvalue of an operator.  Here, the method is probably more
 * accurately called a Picard iteration, since the operator is strictly
 * nonlinear (whereas power iteration implies a linear operator).
 *
 * Acceleration via Steffensen's method is available.  In this case, Aitken's
 * \f$\delta^2\f$ process is used to accelerate \f$k\f$ after two outer
 * iterations.  The Aitken-extrapolants are fed back into response function
 * evaluation, and, in the best case, convergence is second order.
 *
 * The convergence criterion is currently limited to a nonlinear residual 
 * evaluation that makes comparison to Newton methods easier, though other
 * criteria (perhaps on \f$ k \f$ alone) could be used.
 *
 * The outer PowerIter is templated with an \ref InnerIter.  This can be
 * the built-in power iteration (as meant in the traditional sense) or
 * a SLEPc wrapper class allowing power and several Krylov
 * iteration schemes.
 *
 */
//===========================================================================//

template <class Inner>
class PowerIter : public GlobalSolver
{

public:

  /// Typedefs
  //\{
  typedef typename GlobalSolver::SP_globalsolver    SP_base;
  typedef typename util::SP<PowerIter>              SP_globalsolver;
  typedef typename InnerIterBase::SP_inneriter      SP_inneriter;
  //\}

  /*!
   *  \brief Constructor.
   *
   */
  PowerIter(SP_globalproblem problem, SP_globalinput input);

  ///
  ~PowerIter();

  /*!
   *
   * \brief Solve the ERME via power iteration.
   *
   * The inner eigenvalue problem is solved
   *   \f[
   *     \Big ( \mathbf{M}\mathbf{R}(k^{(m)})
   *            - \lambda \mathbf{I} J_{-}^{(n,m)} \Big ) = 0
   *   \f]
   * within an outer iteration \f$ m \f$. The inner iterations can be solved
   * via the power method or one of SLEPc's solvers.
   *
   * The outer iteration is then defined by the eigenvalue update
   *   \f[
   *     k^{(m+1)} = \frac{ \mathbf{F}(k^{(m)})J_{-}^{(n,m)} }
   *                      { \mathbf{L}(k^{(m)})J_{-}^{(n,m)} }
   *   \f]
   * for fission operator \f$ \mathbf{F} \f$ and loss (absorption plus leakage)
   * operator \f$ \mathbf{L} \f$, and where \f$ n \f$ is the number of inner
   * iterations required for convergence on \f$ J_{-} \f$.
   *
   */
  void solve();

private:

  // <-- shouldn't these be private?
  // inherits
  //    GlobalProblem *problem;
  //    vector<scalar> keffHistory;
  //    vector<scalar> keffErrHistory;
  //    vector<scalar> curErrHistory;
  //    vector<scalar> residErrHistory;
  //    scalar wallTime;

  /// \name Data
  //\{

  /// current eigenvalue
  scalar lambda;

  //\}

  /// \name Implementation
  //\{

  /*!
   *  \brief Compute the nonlinear residual for convergence assessment.
   *
   *  \todo Try to combine this and Newton's into one basic resource.
   *
   *  \param    J       Most recent current vector.
   *  \param    keff    Most recent implicit eigenvalue.
   *  \param    lambda  Most recent explicit eigenvalue.
   *  \return           Norm of the nonlinear residual.
   *
   */
  scalar normResid(SP_vector J, scalar keff, scalar lambda);

  // Aitken acceleration
  scalar Aitken(scalar k0, scalar k1, scalar k2);

  //\}

};

#endif // POWERITER_HH
//---------------------------------------------------------------------------//
//                 end of PowerIter.hh
//---------------------------------------------------------------------------//

