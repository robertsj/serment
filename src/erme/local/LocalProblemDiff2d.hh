//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LocalProblemDiff2d.hh
 * \author Jeremy Roberts
 * \date   11/24/2010
 * \brief  Base class for local problem routines.
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 140                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-09-14 12:53:40 -0400 (Wed, 14 Sep 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef LOCALPROBLEMDIFF2D_HH
#define LOCALPROBLEMDIFF2D_HH

#include <string>
#include <vector>
#include "linalg/typedefs.hh"
#include "ResponseFunction.hh"
#include "LocalProblem.hh"
#include "diff2d/Diff2dInput.hh"
#include "diff2d/Diff2dProblem.hh"
#include "diff2d/Diff2dSolver.hh"
#include "diff2d/Diff2dOutput.hh"
using namespace std;

//===========================================================================//
/*!
 * \class LocalProblemDiff2d
 * \brief This class encapsulates the entire Diff2d local problem routine.
 *
 * to be completed
 *
 */
//===========================================================================//

class LocalProblemDiff2d : public LocalProblem
{
  public:
    LocalProblemDiff2d( string in );
    ~LocalProblemDiff2d();
    ResponseFunction **getResponseFunctions( scalar keff );
  private:
    Diff2dInput input;
    vector<Diff2dProblem*> problem;
    vector<Diff2dSolver*> solver;
    ResponseFunction **RFs;
    scalar keffCurrent;
    scalar keffOld;
};

#endif // LOCALPROBLEMDIFF2D_HH

//---------------------------------------------------------------------------//
//                 end of LocalProblemDiff2d.hh
//---------------------------------------------------------------------------//

