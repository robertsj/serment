//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   InvItMatVecWrap.hh
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

#ifndef INVITMATVECWRAP_HH
#define INVITMATVECWRAP_HH
#include "petscmat.h"
#include "LinAlg.hh"
#include "GlobalProblem.hh"

PetscErrorCode InvItMatVecWrap( Mat M, Vec X, Vec Y );

#endif // INVITMATVECWRAP_HH

//---------------------------------------------------------------------------//
//                 end of InvItMatVecWrap.hh
//---------------------------------------------------------------------------//

