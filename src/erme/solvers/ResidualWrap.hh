//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResidualWrap.hh
 * \author Jeremy Roberts
 * \date   10/19/2010
 * \brief  Wrapper for the residual function.
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 140                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-09-14 12:53:40 -0400 (Wed, 14 Sep 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef RESIDUALWRAP_HH
#define RESIDUALWRAP_HH
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
#include "Newton.hh"

PetscErrorCode ResidualWrap( SNES snes, Vec X, Vec F, void *ptr ){
   // cout << " ********** WRAP X ********** " << endl;

   // VecView( X, PETSC_VIEWER_STDOUT_SELF );

    return ((Newton *)ptr)->Residual( snes, X, F, PETSC_NULL );

};

#endif // RESIDUALWRAP_HH

//---------------------------------------------------------------------------//
//                 end of ResidualWrap.hh
//---------------------------------------------------------------------------//

