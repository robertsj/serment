//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   matrix_test.cc
 * \author Jeremy Roberts
 * \date   08/22/2010
 * \brief  test of matrix wrapper for use with petc
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 61                                            $:Rev of last commit
// $Author:: bert                                       $:Author of last commit
// $Date:: 2010-12-05 11:45:26 -0500 (Sun, 05 Dec 2010) $:Date of last commit
//---------------------------------------------------------------------------//

#include <iostream>
#include "petscvec.h"
#include "petscsys.h"
#include "petscksp.h"
#include "SermentVector.hh"
#include "SermentMatrixCRS.hh"
#include "SermentMatrixBCRS.hh"
#include "typedefs.hh"
using namespace std;
//---------------------------------------------------------------------------//

int main(int argc, char *args[])
{
    // initialize petsc
    PetscInitialize(&argc,&args,PETSC_NULL,PETSC_NULL);

	// parameters to make a matrix and two vectors
	integer m, n, nz;
    m  = 10;
    n  = 10;
    nz = 3;

    //-----------------------------------------------------------------------//
    // SermentMatrixCSR Test
    //

  	// create two instances of SermentVector 
    SermentVector X( m ); // we'll fill this one for A*x
	SermentVector Y( m ); // and use this for the result y (i.e. A*x = y)

    // create an instance of SermentMatrixCSR
    SermentMatrixCRS A( m, n, nz ); 

    // create the values and indices to add for the vector
    scalar val;
    val = 1.0;
    scalar vals[10];  
    integer ix[10];
    for (integer i=0; i<X.m; i++)
    {
   		ix[i] = i;
    	vals[i] = val*i;
    }
    X.insertVals( m, ix, vals );
	Y.insertVals( m, ix, vals );

    // create the values and indices to add for the matrix
    scalar      values[3], values2[2];
    integer     numrow, numcol;
    integer     idxrow[1], idxcol[3], idxcol2[2];
    numrow 		= 1;
    values[0] 	= -1.0;  
	values[1] 	= 4.0; 
	values[2] 	= -1.0;

	// build the first and last rows
	integer one = 1.0;integer zero	= 0.0;
    numcol 		= 2; 
    idxrow[0]	= zero; idxcol2[0] 	= zero; idxcol2[1] 	= one;
	values2[0] 	= 2.0; values2[1] 	= -1.0;
    A.insertVals( values2, numrow, idxrow, numcol, idxcol2 );
    idxrow[0]	= m-1; idxcol2[0] 	= m-2; idxcol2[1] 	= m-1;
	values2[0] 	= -1.0; values2[1] 	= 2.0;
    A.insertVals( values2, numrow, idxrow, numcol, idxcol2 );

	// build the rest of the rows
    numcol 		= 3;
    for (integer i=1; i<X.m-1; i++)
    {
        idxrow[0]=i;
        idxcol[0]=i-1; idxcol[1]=i; idxcol[2]=i+1;
        A.insertVals( values, numrow, idxrow, numcol, idxcol );
    }

	// perform A*X = Y
	A.matVec( X, Y );

	// view the matrx
    cout << " A.M = " << endl;
    MatView(A.M,PETSC_VIEWER_STDOUT_SELF);
	// view the vectors
    cout << " X.V = " << endl;
    VecView(X.V,PETSC_VIEWER_STDOUT_SELF);
    cout << " Y.V = " << endl;
    VecView(Y.V,PETSC_VIEWER_STDOUT_SELF);

    // do a linear solve to for verification
    // initialize the solver context
    KSP ksp;
    KSPCreate( PETSC_COMM_WORLD, &ksp );
    KSPSetOperators( ksp, A.M, A.M, DIFFERENT_NONZERO_PATTERN );
    PC prec;
    KSPGetPC(  ksp, &prec );
    KSPSetTolerances( ksp, 1e-5, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT );
    KSPSetFromOptions(ksp);
    KSPSolve( ksp, Y.V, X.V );
    cout << " X.V = " << endl;
    VecView(X.V,PETSC_VIEWER_STDOUT_SELF);
    // which shows that X.V = inv(A.M)Y.V

    //-----------------------------------------------------------------------//
    // SermentMatrixBCSR Test
    //
    // Here, we build in block format the matrix
    //   B = [ C 0 ; 0 D ] 
    // where C is a 3x3 matrix of ones, and D is a 3x3 matrix of twos.
    //
    m = 6;
    n = 1;
    integer bs = 3;
    SermentMatrixBCRS B( m, m, n, bs ); 
    scalar arrayC[9], arrayD[9];
    for (integer i = 0; i<9; i++)
    {
        arrayC[i] = 1.0;
        arrayD[i] = 2.0;
    }
    integer idxr[1];
    integer idxc[1];
    idxr[0] = 0; // add to the 0,0 block element
    idxc[0] = 0;
    B.insertVals( arrayC, n, idxr, n, idxc );
    idxr[0] = 1; // add to the 0,0 block element
    idxc[0] = 1;
    B.insertVals( arrayD, n, idxr, n, idxc );
    B.checkReady();
    cout << " B.M = " << endl;
    MatView(B.M,PETSC_VIEWER_STDOUT_SELF);    

	// destroy the matrix and vectors
    A.releaseMe();
    B.releaseMe();
	X.releaseMe();
	Y.releaseMe();

	PetscFinalize();

	cout << "  ... all finished. now check for errors! " << endl;

    return 0;
}

//---------------------------------------------------------------------------//
//                 end of matrix_test.cc
//---------------------------------------------------------------------------//

