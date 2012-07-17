//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Newton.hh
 * \author Jeremy Roberts
 * \date   11/23/2010
 * \brief  Newton class definition.
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 168                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-11-08 14:23:24 -0500 (Tue, 08 Nov 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef NEWTON_HH
#define NEWTON_HH

#include <iostream>
#include <cmath>
#include "petscsnes.h"
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
#include "JacobianShell.hh"

//===========================================================================//
/*!
 * \class Newton
 * \brief This class solves GlobalProblems via use of Newton's method.
 *
 * to be completed
 *
 */
//===========================================================================//

class Newton : public GlobalSolver
{

public:

  /// Typedefs
  //\{
  typedef typename GlobalSolver::SP_globalsolver SP_base;
  typedef typename util::SP<Newton> SP_globalsolver;
  typedef typename GlobalProblem::SP_globalproblem SP_globalproblem;
  typedef typename GlobalInput::SP_globalinput SP_globalinput;
  //\}

  Newton(SP_globalproblem problem, SP_globalinput input);

  ~Newton();

  void solve();

  PetscErrorCode Residual(SNES snes, Vec X, Vec F, void *ptr);

  //protected:
  // inherits
  //    GlobalProblem *problem;
  //    vector<scalar> keffHistory;
  //    vector<scalar> keffErrHistory;
  //    vector<scalar> curErrHistory;
  //    vector<scalar> residErrHistory;
  //    scalar wallTime;
  SermentMatrix *P; // the preconditioner
  JacobianShell *fp; // the jacobian
  SermentVector *x; // the unknown vector
  SermentVector *f; // the nonlinear residual vector
  SNES snes; /* SNES context */

  scalar lambda;// "current eigenvalue"
};

#endif // NEWTON_HH
//---------------------------------------------------------------------------//
//                 end of Newton.hh
//---------------------------------------------------------------------------//

