//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SermentMatrixFull.cc
 * \author Jeremy Roberts
 * \date   10/22/2010
 * \brief  Member definitions of abstract class SermentMatrixFull
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 68                                            $:Rev of last commit
// $Author:: bert                                       $:Author of last commit
// $Date:: 2010-12-09 17:17:34 -0500 (Thu, 09 Dec 2010) $:Date of last commit
//---------------------------------------------------------------------------//
#include <iostream>
#include "petscvec.h"
#include "petscmat.h"
#include "SermentMatrixFull.hh"
#include "SermentVector.hh"
#include "typedefs.hh"


//---------------------------------------------------------------------------//
/*!
 * \brief Matrix-vector multiplication, Mx-->y 
 *
 * This function performs matrix-vector multiplcication of the form y=Mx.  The
 * actual multiplication is done by the Petsc MatMult function, where x and y
 * cannot be the same.
 */
void SermentMatrixFull::matVec( SermentVector &x, SermentVector &y )
{
    x.checkReady();
	y.checkReady();
	this->checkReady();
    MatMult( M, x.V, y.V );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Matrix-vector multiplication, Mx-->x 
 *
 * This function performs matrix-vector multiplcication of the form x=Mx.  The
 * actual multiplication is done by the Petsc MatMult function.  Because Petsc's
 * function does not allow x and y to be the same, we initiate a local Petsc
 * with vector into which the results are placed.  The contents of x.V are then
 * replaced with the output.
 */
void SermentMatrixFull::matVec( SermentVector x )
{
  //  Vec     y;          // temporary Petsc vector
  //  scalar  xvals;      // scalar array for values of x.V
    
   // MatMult( M, x.V, y )
    // define me later
    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Method check whether SermentVector is ready for operations
 *
 * Because the underlying external library (here Petsc) has special operations
 * to ready vectors for operations, they are inserted here as a check to be
 * called before any operation.  If isReady is true, nothing happens, but if
 * it is false, the vector is assembled.
 */
void SermentMatrixFull::checkReady()
{
	if ( isReady == false ) 
	{
		MatAssemblyBegin( M, MAT_FINAL_ASSEMBLY );
		MatAssemblyEnd( M, MAT_FINAL_ASSEMBLY );
        isReady = true;
	}
	return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief This function displays the matrix contents.
 *
 */
void SermentMatrixFull::viewMe()
{
	checkReady();
    MatView( M, PETSC_VIEWER_STDOUT_SELF );
	return;
}


//---------------------------------------------------------------------------//
//                 end of SermentMatrixFull.cc
//---------------------------------------------------------------------------//

