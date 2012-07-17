//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Diff2dSolver.hh
 * \author Jeremy Roberts
 * \date   10/26/2010
 * \brief  A class for solving various diff2d problems.
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 171                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-11-19 20:36:52 -0500 (Sat, 19 Nov 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef DIFF2DSOLVER_HH
#define DIFF2DSOLVER_HH

#include <iostream>
#include <fstream>
#include <vector>
#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"
#include "petscpc.h"
#include "petscsys.h"
#include "../linalg/typedefs.hh"
#include "Diff2dInput.hh"
#include "Diff2dProblem.hh"
#include "utilities/SP.hh"
using namespace std;

//===========================================================================//
/*!
 * \class Diff2dSolver
 * \brief A class for solving diff2d problems.
 *
 *  Diff2dSolver is blah blah blah.
 *
 */
//===========================================================================//

class Diff2dSolver
{
public:

  typedef util::SP<Diff2dSolver> SP_solver;

  // the flux
  vector<Vec> phi;
  Diff2dSolver(){}
  Diff2dSolver(Diff2dInput &input, Diff2dProblem &problem, int el);
  ~Diff2dSolver();
  void solve();
  void compRespFct();
  void destroy();
private:
  Diff2dInput *in;
  Diff2dProblem *pr;
  integer elid;
  KSP ksp;
  PC prec;
  Vec scsrc, fisrc, tmp, fisrc0;
  vector<Vec> phi0;
  integer maxit, nrow, ptype;
  scalar epsk, epss, keff;
  void solveFixedNoMult()
  {
    return;
  }
  ;
  void solveFixedMult();
  void solveEigenvalue();
};

#endif // DIFF2DSOLVER_HH
//---------------------------------------------------------------------------//
//                 end of Diff2dSolver.hh
//---------------------------------------------------------------------------//

