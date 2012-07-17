//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   JacobianEmpty.hh
 * \author Jeremy Roberts
 * \date   10/19/2010
 * \brief  Wrapper for the Jacobian's MyMatMult function.
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 140                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-09-14 12:53:40 -0400 (Wed, 14 Sep 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef JACOBIANEMPTY_HH
#define JACOBIANEMPTY_HH
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
#include "JacobianShell.hh"

PetscErrorCode JacobianEmpty( SNES snes, Vec X, Mat *J, Mat *B, MatStructure *flag, void *ptr  )
{
   // Newton *me = (Newton*)ptr;
    //VecCopy( X, me->x->V );
    return 0;
};

#endif // JACOBIANEMPTY_HH

//---------------------------------------------------------------------------//
//                 end of JacobianEmpty.hh
//---------------------------------------------------------------------------//

