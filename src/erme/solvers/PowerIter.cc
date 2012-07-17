//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PowerIter.cc
 * \author Jeremy Roberts
 * \date   11/24/2010
 * \brief  Member definitions of base class PowerIter
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 177                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-12-09 16:31:49 -0500 (Fri, 09 Dec 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#include <iostream>
#include <cmath>

#include "serment_config.h"

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
#include "PowerIter.hh"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
template <class Inner>
PowerIter<Inner>::PowerIter(SP_globalproblem problem, SP_globalinput input) :
  GlobalSolver(problem, input)
{
  std::cout << " CONSTRUCTING PowerIter " << std::endl;
  // nothing more here for now
}

//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
template <class Inner>
PowerIter<Inner>::~PowerIter()
{
  // nothing more here right now, as M is destroyed explicitly in Release()
  return;
}

// Solve the ERME via power iteration.
template <class Inner>
void PowerIter<Inner>::solve()
{

  // User-provided initial guess
  scalar keff = d_input->keff;

  // Hard-coded initial guess
  double lambda = 1.0;

  // Temporary working vectors.
  SP_vector J0;
  SP_vector J1;
  J0 = new SermentVector(d_input->degfree);
  J1 = new SermentVector(d_input->degfree);
  // and dummy place holder.
  SP_vector Jtmp;

//  if (!d_problem->run_already()) \todo Give the solver the final solution
//  {
    // Uniform zeroth order initial guess
    J0->vecSet(0.0);
    for (integer i = 0; i < d_input->degfree; i = i + (d_input->spaceord + 1)
        * (d_input->angleord + 1) * d_input->numgroups)
    {
      J0->insertVal(i, 1.0);
    }
    // Normalize it.
    lambda = sqrt(J0->vecDot(*J0));
    J0->vecScale(1.0 / lambda);
//  }
//  else
//  {
//    // Use the last solution.
//    d_problem->J.vecCopy(*J0);
//  }

  scalar loss       = 0;
  scalar gain       = 0;
  scalar absorption = 0;
  scalar leakage    = 0;
  scalar keffo      = 0;
  scalar keffoo     = 0;
  scalar lambdao    = 0;
  scalar lambdaoo   = 0;

  bool aitken       = false;     // Used only to view Aitken estimate of k
  bool steffensen   = false;     // Actually plugs Aitken estimate in
  if (d_input->ilulevel == 1)
    steffensen = true;
  else
    aitken = true;

  integer innertot  = 0; // count total iterations
  scalar  norm      = 0; // residual norm


  cout << " beginning power iterations.... ";
  if (steffensen) cout << " with Steffensen ";
  cout << endl;

  // Initial residual norm
  norm = normResid(J0, keff, lambda);
  printf(
      "  %3i  POWER ITERATION function norm  %8.6e lambda %12.9f keff %12.9f keffAit %12.9f cumulative inners %8i \n",
      0, norm, 1.0, keff, 1.0, 0);

  // Create inner solver.
  Inner solver(d_problem->M, d_problem->R);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
  // outers
  int it; // count outer iteration
  for (it = 1; it <= d_input->maxit; it++)
  {

    // Ensure all response operators are up-to-date.
    d_problem->R->updateData(keff);
    d_problem->L.updateData(keff);
    d_problem->A.updateData(keff);
    d_problem->F.updateData(keff);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
    // inners

    lambda = solver.solve(d_input->maxit, d_input->epss, J0, J1);

    // end inners
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

    // Compute gains and losses.
    gain        = d_problem->F.vecDot(*J1);
    absorption  = d_problem->A.vecDot(*J1);
    leakage     = d_problem->L.computeLeakage(*J1);
    loss        = absorption + leakage;

    // Keep oldest and second oldest keff values
    keffoo  = keffo;
    keffo   = keff;
    // and update keff.
    keff    = gain / loss;
    // Do the same for lambda.
    lambdaoo  = lambdao;
    lambdao   = lambda;
    // If requested, use Aitken extrapolant, yielding Steffensen's method.
    if (it > 2 and steffensen)
    {
      keff      = Aitken(keffoo, keffo, keff);
      lambda    = Aitken(lambdaoo, lambdao, lambda);
    }


    // Compute the norm of the nonlinear residual.
    norm = normResid(J1, keff, lambda);

    // Save current vector.
    Jtmp = J1;
    J1   = J0;
    J0   = Jtmp;

    printf("  %3i  POWER ITERATION function norm  %8.6e lambda %12.9f keff %12.9f keffAit %12.9f cumulative inners %8i \n",
           it, norm, lambda, keff, Aitken(keffoo, keffo, keff), innertot);

    // test fissionrates
    //fissionRates(J0);

    // the norm of the residual is sufficiently small
    if (norm < d_input->epss) break;

    // warning about convergence
    if (it == d_input->maxit - 1) cout
        << "****warning, reached max outers w/o convergence" << endl;

  }
  // end outers loop
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

  std::printf(" FINAL POWER ITERATION EIGENVALUES: \n");
  std::printf(" **** FINAL KEFF        = %12.9f \n", keff);
  std::printf(" **** FINAL LAMBDA      = %12.9f \n", lambda);
  std::printf(" **** OUTER ITERATIONS  = %8i \n", it);
  std::printf(" **** INNER ITERATIONS  = %8i \n", innertot);
 // std::printf(" **** LINEAR ITERATIONS = %8i \n", linits);

  // Copy J0 to final solution
  //J0->vecCopy(d_problem->J);
  //d_problem->set_run_already();
  //  J0->viewMe();

  return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief This function computes the \f$ L_2 \f$ norm of the nonlinear residual.
 *
 * The nonlinear residual is defined as
 *   \f[
 *       \mathbf{f(x)} = \left [\begin{array}{c}
 *	        (\mathbf{M}\mathbf{R}(k)-\lambda \mathbf{I}) \mathbf{J_-} \\
 *	        \mathbf{F}(k)\mathbf{J_-} - (k\mathbf{L}(k)\mathbf{J_-} ) \\
 *	        \frac{1}{2} \mathbf{J^T_-} \mathbf{J_-} - \frac{1}{2}  
 *	      \end{array} 
 *       \right ]  = \mathbf{0} \, ,
 *   \f]
 * which is the same as used in the Newton-based schemes.  The \f$ L_2 \f$
 * norm is then \f$ \sqrt{ \mathbf{f(x)}^T \mathbf{f(x)} } \f$.
 *
 */
template <class Inner>
scalar PowerIter<Inner>::normResid(SP_vector J,
                                   scalar    keff,
                                   scalar    lambda)
{

  // update response quantities
  d_problem->R->updateData(keff);
  d_problem->L.updateData(keff);
  d_problem->A.updateData(keff);
  d_problem->F.updateData(keff);

  // Fj, Fk, and Flambda are the residuals
  SermentVector Fj(d_input->degfree);
  SermentVector temp(d_input->degfree);

  d_problem->R->matVec(*J, temp);
  d_problem->M->matVec(temp, Fj);

  Fj.vecAYPV(-lambda, *J); // = (M*R-lambda*I)J

  scalar normFjsq = Fj.vecDot(Fj); // = ||(M*R-lambda*I)J||^2
  scalar normFk = d_problem->F.vecDot(*J) - // = (F-kL)J
      keff * (d_problem->A.vecDot(*J) + d_problem->L.computeLeakage(*J));
  scalar normFlambda = 0.5 - 0.5 * J->vecDot(*J);
  scalar norm = sqrt(normFjsq + pow(normFk, 2.0) + pow(normFlambda, 2.0));

#ifdef DEBUG
  cout << " normFj      = " << sqrt( normFjsq ) << endl;
  cout << " normFk      = " << normFk << endl;
  cout << " normFlambda = " << normFlambda << endl;
#endif

  Fj.releaseMe();
  temp.releaseMe();

  return norm;
}

//---------------------------------------------------------------------------//
/*!
 * \brief This function computes the Aitken/Steffensen estimate.
 *
 * Aitken's \f$ \Delta^2 \f$ method is an extremely easy way to
 * accelerate convergence of a series satisfying a few constraints.  Suppose we
 * have a sequence of estimates for the eigenvalue: 
 * \f$ k^{(n-2)} \f$, \f$ k^{(n-1)} \f$, and \f$ k^{(n)} \f$, which converge to 
 * \f$ k^* \f$ as \f$ n\to \infty \f$.  We define a new sequence
 * \f$ k'^{(n)} \f$ such that
 * \f[
 k'^{(n)} = k^{(n)} - \frac{ ( k^{(n)} - k^{(n-1)} )^2 }{ k^{(n)} - 2k^{(n-1)} + k^{(n-2)} } \, .
 * \f]
 * It can be shown that this new sequence \f$ k'^{(n)} \f$  converges to 
 * \f$ k^* \f$ faster than the sequence \f$ k^{(n)} \f$ for monotonically 
 * converging \f$ |k^{(n)}| \f$.  Moreover, this new update can be used to 
 * update the response quantities; used in this way, Aitken's becomes 
 * Steffensen's method.
 *
 */
template <class Inner>
scalar PowerIter<Inner>::Aitken(scalar k0, scalar k1, scalar k2)
{
  scalar kA = k0 - pow(k1 - k0, 2) / (k2 - 2.0 * k1 + k0);
  return kA;
}

// Explicit instantiations;
template class PowerIter<InnerIterPower> ;
#ifdef SERMENT_ENABLE_SLEPC
template class PowerIter<InnerIterSLEPc> ;
#endif
//---------------------------------------------------------------------------//
//                 end of PowerIter.cc
//---------------------------------------------------------------------------//

