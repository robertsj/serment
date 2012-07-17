//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   vector_test.cc
 * \author Jeremy Roberts
 * \date   08/10/2010
 * \brief  test of vector wrapper for use with petc
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 20                                            $:Rev of last commit
// $Author:: bert                                       $:Author of last commit
// $Date:: 2010-10-26 10:23:31 -0400 (Tue, 26 Oct 2010) $:Date of last commit
//---------------------------------------------------------------------------//

#include <iostream>
#include "petscvec.h"
#include "petscsys.h"
#include "SermentVector.hh"
#include "SermentMatrixFull.hh"
#include "typedefs.hh"
using namespace std;
//---------------------------------------------------------------------------//

int main(int argc, char *args[])
{
    // initialize petsc
    PetscInitialize(&argc,&args,PETSC_NULL,PETSC_NULL);

	// parameters to make a matrix of size...
	integer m;
    m = 10;

  	// create an instance of SermentVector 
    SermentVector A(m);

    // now insert values using the single value approach
    scalar val;
    val = 1.0;
	for (integer i=0; i<A.m; i++) 
    {
		A.insertVal( i, val*i );
	}

	// finalize the vector after construction
	VecAssemblyBegin(A.V);
	VecAssemblyEnd(A.V);

	// view the matrix
    cout << " A.V = " << endl;
    VecView(A.V,PETSC_VIEWER_STDOUT_SELF);

    // now do the same for a new vector using InsertVals
    SermentVector B(m);
    
    // create the values and indices to add
    scalar vals[10];
    integer ix[10];
    for (integer i=0; i<B.m; i++)
    {
   		ix[i] = i;
    	vals[i] = val*i;
    }
    B.insertVals( m, ix, vals );

	// finalize the matrix after construction
	VecAssemblyBegin(B.V);
	VecAssemblyEnd(B.V);

	// view the matrix
    cout << " B.V = " << endl;
    VecView(B.V,PETSC_VIEWER_STDOUT_SELF);


    return 0;

	PetscFinalize();
}

//---------------------------------------------------------------------------//
//                 end of vector_test.cc
//---------------------------------------------------------------------------//

