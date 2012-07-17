//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SermentVector.cc
 * \author Jeremy Roberts
 * \date   10/15/2010
 * \brief  Member definitions of class SermentVector
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 177                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-12-09 16:31:49 -0500 (Fri, 09 Dec 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#include <iostream>
#include "petscvec.h"
#include "SermentVector.hh"
#include "typedefs.hh"
#include "../utilities/GenException.hh"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
SermentVector::SermentVector(integer a)
{
	m = a;
    isReady = false;
	VecCreateSeq( PETSC_COMM_SELF, m, &V );
}

//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
SermentVector::~SermentVector()
{
    return; // putting destruction of Petsc Vec in its own routine
}


//---------------------------------------------------------------------------//
/*!
 * \brief Method to insert single indexed value into the vector. 
 *
 * Note, this is <b>not</b> the ideal way to construct Petsc vectors and is
 * included primarily for testing.
 */
void SermentVector::insertVal( integer row, scalar value )
{
	VecSetValue( V, row, value, INSERT_VALUES );
    isReady = false; // change values = no longer ready
	return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Method to insert several indexed values into the matrix. 
 *
 * For standard compressed row storage matrices, this is the preferred method
 * for inserting values.  Arguments:
 *   ni - number of elements to add
 *   ix - indices where to add
 *    y - values to add
 * The Petsc documentation is as follows:
 *	x 	- vector to insert in
 *	ni 	- number of elements to add
 *	ix 	- indices where to add
 *	y 	- array of values
 *	iora 	- either INSERT_VALUES or ADD_VALUES, where ADD_VALUES adds values 
 *  to any existing entries, and INSERT_VALUES replaces existing entries with 
 *  new values 
 */
void SermentVector::insertVals( integer ni, const integer ix[],  
                                const scalar y[] )
{
    VecSetValues( V, ni, ix, y, INSERT_VALUES );
	isReady = false; // make sure after insertions that construction is done
	return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set all entries to single value.
 *
 */
void SermentVector::vecSet( scalar a )
{
    VecSet(V, a);
    isReady = false; // change values = no longer ready
}

//---------------------------------------------------------------------------//
/*!
 * \brief V = a*V + Y 
 */
void SermentVector::vecAVPY( scalar a, SermentVector &Y )
{
	Y.checkReady();
	this->checkReady();
	// y = alpha*y + x
    VecAYPX( V, a, Y.V );
	return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief V = a*Y + V 
 */
void SermentVector::vecAYPV( scalar a, SermentVector &Y )
{
	Y.checkReady();
	this->checkReady();
	// y = alpha*x + y
    VecAXPY( V, a, Y.V);
	return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief val = y'x  
 */
scalar SermentVector::vecDot( SermentVector &Y )
{
	Y.checkReady();
	this->checkReady();
	// VecDot(Vec x,Vec y,PetscScalar *val)
    scalar val;
    VecDot( V, Y.V, &val );
	return val;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Z = X.*Y
 */
void SermentVector::vecPointMult( SermentVector &Y, SermentVector &Z )
{
	Y.checkReady();
	Z.checkReady();
	this->checkReady();

	// temporary arrays
	scalar *x_a, *y_a, *z_a;

	VecGetArray( this->V, &x_a );
	VecGetArray( Y.V,     &y_a );
	VecGetArray( Z.V,     &z_a );    

	// ensure they are the same size
	if ( this->Length() != Z.Length() or this->Length() != Y.Length() )
	{
		std::cout << __LINE__ << " in " << __FILE__ << std::endl;
		throw GenException(__LINE__,__FILE__,"Inappropriate vector lengths.");
	}

	// loop through and compute
    for ( int i = 0; i < this->Length(); i++ )
	{
		z_a[i] = x_a[i] * y_a[i];
	}

 	VecRestoreArray( this->V, &x_a );
 	VecRestoreArray( Y.V,     &y_a );
 	VecRestoreArray( Z.V,     &z_a );

	return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief x = a*x
 *
 * \param a   scalar value by which the vector x is scaled
 */
void SermentVector::vecScale( scalar a )
{
	this->checkReady();
	// VecScale (Vec x, PetscScalar alpha)
    VecScale (V, a);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Copy me to Y
 *
 * \param a   scalar value by which the vector x is scaled
 */
void SermentVector::vecCopy( SermentVector &Y )
{
    this->checkReady();
    // VecCopy(Vec x,Vec y)
    VecCopy (V, Y.V);
}


//---------------------------------------------------------------------------//
/*!
 * \brief Return my length.
 */
integer SermentVector::Length()
{
	return this->m;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Method to release memory for the Petsc Vec
 *
 * Ideally, this method would be part of the destructor, but since it would
 * be called after PetscFinalize---which seems to do deallocations of its
 * own---it's best to use explicit deallocation when a SermentVector is no
 * longer needed.
 */
void SermentVector::releaseMe()
{
	VecDestroy( &V );
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
void SermentVector::checkReady()
{
	if ( isReady == false ) 
	{
		VecAssemblyBegin( V );
		VecAssemblyEnd( V );
        isReady = true;
	}
	return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief View the vector contents. 
 *
 */
void SermentVector::viewMe()
{
    checkReady();
	VecView( V,PETSC_VIEWER_STDOUT_SELF );
	return;
}


//---------------------------------------------------------------------------//
//                 end of SermentVector.cc
//---------------------------------------------------------------------------//

