//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   InvItShell.cc
 * \author Jeremy Roberts
 * \date   10/25/2010
 * \brief  Member definitions of abstract class InvItShell
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 168                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-11-08 14:23:24 -0500 (Tue, 08 Nov 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#include <iostream>
#include "LinAlg.hh"
#include "GlobalProblem.hh"
#include "InvItMatVecWrap.hh"
#include "InvItShell.hh"

//---------------------------------------------------------------------------//
/*!
 * \brief This performs the action of (M*R-lambda*I) for use in inverse it.
 *
 */
void InvItShell::myMatVec( Mat &M, Vec &x, Vec &y )
{

//    problem->R.updateData( k ); // make sure up-to-date
//    SermentVector temp( m );
//
//    MatMult( problem->R.M, x, temp.V );  // = R*x
//    MatMult( problem->M.M, temp.V, y );  // = M*R*x
//    VecAXPY( y, -lambda, x );            // = M*R*x - lambda*x
//
//    temp.releaseMe();
//
//    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief This is a matrix-vector multiplication wrapper function.
 *
 * It exists outside the SermentMatrixShell class so that it can be passed
 * to Petsc, as Petsc needs pointers to functions, and it doesn't like 
 * pointers to member methods.
 *
 */
PetscErrorCode InvItMatVecWrap( Mat M, Vec f, Vec fpf )
{

    void *ctx;
    MatShellGetContext( M, &ctx );
    InvItShell *me = (InvItShell*)ctx;

    me->myMatVec( M, f, fpf );

    return 0;
}

//---------------------------------------------------------------------------//
//                 end of InvItShell.cc
//---------------------------------------------------------------------------//

