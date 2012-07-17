//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   InnerIterPower.cc
 * \author Jeremy Roberts
 * \date   Nov 7, 2011
 * \brief  InnerIterPower member definitions.
 * \note   Copyright (C) 2011 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//
// $Rev::                                               $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date::                                              $:Date of last commit
//---------------------------------------------------------------------------//

#include <cmath>
#include <iostream>

#include "InnerIterPower.hh"

#include "utilities/DBC.hh"

#include "serment_config.h"

InnerIterPower::InnerIterPower(SP_M M, SP_R R)
  : InnerIterBase(M, R)
{
  // Nothing here for now.
}

// Solve M*R*J = lambda*J for the dominant eigenpair.
scalar InnerIterPower::solve(int            max_iters,
                             double         tol,
                             SP_vector      J_in,
                             SP_vector      J)
{
  Require(nax_iters > 0);
  Require(tol > 0);
  Require(J_in);            // No null
  Require(J);               // vectors

  // Set up temporary pointers to cycle around the allocated storage.
  //   Remember, the client has storage in J_in and J, and we have
  //   storage in our own tmp vector.  The third vector is needed
  //   because two matrices are applied and we DO NOT change the
  //   incoming vector.
  SP_vector J1, J2, J3, Jtmp;
  J1 = J_in;
  J2 = d_J_tmp;
  J3 = J;

  // Eigenvalue and error.
  scalar lambda_0   = 0.0;
  scalar lambda     = 0.0;
  scalar lambda_err = 0.0;

  // Eigenvector errors.  We'll estimate the dominance ratio using
  //   a few iterations, and then estimate how many iterations until
  //   the error has died out such that the tolerance is met.  We
  //   use J_err = ||MRJ-lambda*J||, as this is the bulk of the
  //   nonlinear residual; get this right, get f(x) right.
  scalar J_err_0  = 1.0;
  scalar J_err    = 1.0;

  // Make sure J_in is normalized to one.
  lambda = std::sqrt(J1->vecDot(*J1));
  J1->vecScale(1.0 / lambda);

  // Next iteration to check eigenvector convergence.
  int check_err = true;

  for (int i = 0; i < max_iters; i++)
  {

    // Save previous value.
    lambda_0 = lambda;
    J_err_0  = J_err;

    // Apply the MR operator.
    d_R->matVec(*J1, *J2);
    d_M->matVec(*J2, *J3);

    // Compute the current eigenvalue and update eigenvalue error.
    lambda     = sqrt(J3->vecDot(*J3));
    lambda_err = std::fabs(lambda - lambda_0);

    // If needed, update eigenvector residual and compute next check point.
    if (i < 5 or i % 200 == 0)
    {

      // J3 = M*R*J1, where ||J1|| = 1. Hence, the res is J3-J1.
      J1->vecAVPY(-lambda, *J3);        // J1   <-- J1 - J3
      J_err = sqrt(J1->vecDot(*J1));    // Jerr <-- ||J1 - J3|| = norm(res)
      scalar rho   = J_err / J_err_0;   // rho ~ ||res(n)|| / ||res(-1)||
      Ensure(rho < 1.0)                 // Stops divergence.

      // Using the fact J_err(N) ~ J_err(0)*rho^N, we check again in
      //   N = ln(tol/J_err)-ln(rho) iterations.
      check_err = min(300, int(std::log(tol/J_err) / std::log(rho)));
      Ensure(check_err >= 0);
      check_err += i;
#ifdef DEBUG
      std::cout << " Dominance ratio estimated to be: " << rho
                << ". Checking eigenvector error next at i = " << check_err
                << ". Error is " << J_err << " tol is " << tol << std::endl;
#endif

    }

    // Normalize the current.
    J3->vecScale(1.0 / lambda);

    // Cycle the pointers.
    Jtmp = J3;
    J3   = J2;
    J2   = J1;
    J1   = Jtmp;

    if ( i % 200 == 0 )
         std::printf("  i %3i  lambda %12.9f  lambda_err %12.9e \n",
                     i, lambda, lambda_err);

    // Check for convergence.
    if (J_err < tol and lambda_err < tol)
    {
      // Copy most recent eigenvector to the intended output memory.  If
      //   J1 and J actually point to the same memory already, the PETSc
      //   copy function knows that and simply returns.
      J1->vecCopy(*J);
      // Better have a postive eigenvalue!  Otherwise, an error lurks.
      Ensure(lambda >= 0.0);
      std::printf("Successful inners convergence in %5i iterations \n", i);
      return lambda;
    }
  }
  std::printf("Warning: inners failed to converge.  J error = %12.9f \n",
              J_err);
  return lambda;
}

//---------------------------------------------------------------------------//
//              end of InnerIterPower.cc
//---------------------------------------------------------------------------//
