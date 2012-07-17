//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SermentMatrixShell.cc
 * \author Jeremy Roberts
 * \date   10/25/2010
 * \brief  Member definitions of abstract class SermentMatrixShell
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 70                                            $:Rev of last commit
// $Author:: bert                                       $:Author of last commit
// $Date:: 2010-12-10 18:11:28 -0500 (Fri, 10 Dec 2010) $:Date of last commit
//---------------------------------------------------------------------------//

#include <iostream>
#include "petscvec.h"
#include "petscmat.h"
#include "SermentMatrixShell.hh"
#include "SermentVector.hh"
#include "typedefs.hh"

using namespace std;

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
SermentMatrixShell::SermentMatrixShell( integer a, integer b, void *ctx ) 
      : SermentMatrix(a,b)
{
	MatCreateShell( PETSC_COMM_SELF, m, n, PETSC_DETERMINE, PETSC_DETERMINE, 
                    ctx, &M );
}


//---------------------------------------------------------------------------//
/*!
 * \brief Matrix-vector multiplication, Mx-->y 
 *
 * This function performs matrix-vector multiplcication of the form y=Mx.  The
 * actual multiplication is done by the MyMatVec function, allowing for a 
 * completely matrix-free approach.
 */
void SermentMatrixShell::matVec( SermentVector &x, SermentVector &y )
{
    x.checkReady();
	y.checkReady();

    myMatVec( M, x.V, y.V );

}


//---------------------------------------------------------------------------//
//                 end of SermentMatrixShell.cc
//---------------------------------------------------------------------------//

