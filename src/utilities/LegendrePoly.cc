//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LegendrePoly.cc
 * \author Jeremy Roberts
 * \date   10/26/2010
 * \brief  Member definitions of class LegendrePoly
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 177                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-12-09 16:31:49 -0500 (Fri, 09 Dec 2011) $:Date of last commit
//---------------------------------------------------------------------------//
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "petscvec.h"
#include "LegendrePoly.hh"
#include "../linalg/typedefs.hh"

using namespace std;

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
LegendrePoly::LegendrePoly()
{
    // no constructor behavior
}

//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
LegendrePoly::~LegendrePoly()
{
    // nothing for destruction
}

//---------------------------------------------------------------------------//
/*!
 * \brief This deallocates the Petsc Vecs explicitly.
 */
void LegendrePoly::destroy()
{
    for (integer i = 0; i <= maxOrder; i++)
    {
        VecDestroy( &Px[i] );
        VecDestroy( &Py[i] );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief This creates the discrete Legendre polynomials for all edges.
 *
 */
void LegendrePoly::buildMe( integer nx, integer ny, integer mxord )
{
    integer Nx = nx - 1; // max degree for horizontal surface
    integer Ny = ny - 1; // max degree for vertical surface
    numx = nx; numy = ny; maxOrder = mxord;
    if ( maxOrder+1 > nx )
        cout << " WARNING: LP maxOrder exceeds X fine mesh count!" << endl;
    if ( maxOrder+1 > ny )
        cout << " WARNING: LP maxOrder exceeds Y fine mesh count!" << endl;

    integer kx[nx], ky[ny];
    for ( integer i = 0; i < nx; i++ )
        kx[i] = i;
    for ( integer i = 0; i < ny; i++ )
        ky[i] = i;

    scalar one = 1.0;

    // initialize the Vecs
    Px.resize(maxOrder+1);
    Py.resize(maxOrder+1);
    VecCreate( PETSC_COMM_SELF, &Py[0] );
    VecCreate( PETSC_COMM_SELF, &Px[0] );
    VecSetSizes( Px[0], PETSC_DECIDE, nx );
    VecSetSizes( Py[0], PETSC_DECIDE, ny );
    VecSetFromOptions( Px[0] );
    VecSetFromOptions( Py[0] );
    VecSet( Px[0], one );
    VecSet( Py[0], one );
    VecAssemblyBegin( Px[0] ); VecAssemblyEnd( Px[0] );
    VecAssemblyBegin( Py[0] ); VecAssemblyEnd( Py[0] );

    for ( integer i = 1; i < maxOrder+1; i++ )
    {
        VecDuplicate( Px[0], &Px[i] );
        VecDuplicate( Py[0], &Py[i] );    
    }

    // now build the actual vectors using arrays and then VecSetValues
    // where kx and ky are used as the indices
    scalar tmpPxA[nx], tmpPxB[nx], tmpPxC[nx];
    scalar tmpPyA[ny], tmpPyB[ny], tmpPyC[ny];
    scalar *previous, *current, *newest, *temporary;

    // P0 is set; now generate P0 and P1 in temporary array, set P1, and then
    // proceed to build P3, P4, ...
    for ( integer i = 0; i < nx; i++ )
    {
        tmpPxA[i] = one;
        tmpPxB[i] = 1.0 - 2.0*i/Nx;
    }
    if ( maxOrder > 0 )
    {
        VecSetValues( Px[1], nx, kx, tmpPxB, INSERT_VALUES );
        VecAssemblyBegin( Px[1] ); 
        VecAssemblyEnd( Px[1] );
    }
    previous = &tmpPxA[0]; current = &tmpPxB[0]; newest = &tmpPxC[0];
    scalar ii;
    for ( integer i = 2; i < maxOrder+1; i++ )
    {
        ii=i-1.0;
        for ( integer j = 0; j < nx; j++ )
        {
            newest[j] = ( (2.0*ii+1)*(Nx-2.0*j) * current[j] - 
                        ii*(Nx+ii+1.0) * previous[j] ) / ((ii+1.0)*(Nx-ii)) ;
        }
        VecSetValues( Px[i], nx, kx, newest, INSERT_VALUES );
        VecAssemblyBegin( Px[i] ); 
        VecAssemblyEnd( Px[i] );
        temporary   = newest; 
        newest      = previous; 
        previous    = current; 
        current     = temporary;
    }

    // P0 is set; now generate P0 and P1 in temporary array, set P1, and then
    // proceed to build P3, P4, ...
    for ( integer i = 0; i < ny; i++ )
    {
        tmpPyA[i] = one;
        tmpPyB[i] = 1.0 - 2.0*i/Ny;
    }
    if ( maxOrder > 0 )
    {
        VecSetValues( Py[1], ny, ky, tmpPyB, INSERT_VALUES );
        VecAssemblyBegin( Py[1] ); 
        VecAssemblyEnd( Py[1] );
    }

    previous = &tmpPyA[0]; current = &tmpPyB[0]; newest = &tmpPyC[0];

    for ( integer i = 2; i < maxOrder+1; i++ )
    {
        ii=i-1.0;
        for ( integer j = 0; j < ny; j++ )
        {
            newest[j] = ( (2.0*ii+1)*(Ny-2.0*j) * current[j] - 
                        ii*(Ny+ii+1.0) * previous[j] ) / ((ii+1.0)*(Ny-ii)) ;
        }
        VecSetValues( Py[i], ny, ky, newest, INSERT_VALUES );
        VecAssemblyBegin( Py[i] ); 
        VecAssemblyEnd( Py[i] );
        temporary = newest; 
        newest = previous; 
        previous = current; 
        current = temporary;
    }

    makeWeights();
}

//---------------------------------------------------------------------------//
/*!
 * \brief This function creates the normalization weights.
 *
 */
void LegendrePoly::makeWeights()
{
    wx.resize(maxOrder+1);
    wy.resize(maxOrder+1);
    // do for x
    scalar a,b,c;
    for ( integer i = 0; i <= maxOrder; i++ )
    {
        a = 2*i+1;
        b = factorial(numx+i)/factorial(numx-1);
        c = factorial(numx-1)/factorial(numx-1-i);
        wx[i] = a*c/b;
    }
    // do for y
    for ( integer i = 0; i <= maxOrder; i++ )
    {
        a = 2*i+1;
        b = factorial(numy+i)/factorial(numy-1);
        c = factorial(numy-1)/factorial(numy-1-i);
        wy[i] = a*c/b;
    }

}

//---------------------------------------------------------------------------//
/*!
 * \brief This function expands a partial current and gives back coefficients.
 *
 * The user gives the partial current Petsc Vec and an array to be filled with
 * the coefficients.  Note, the array *must* be the correct size, i.e. order+1,
 * as there is currently no check.  Example use: scalar coefsX[order+1], 
 * LP.expandCur( PcurX, coefsX );
 */
void LegendrePoly::expandCur( Vec &partialCurrent, scalar legCoefs[] )
{
    integer size;
    scalar tmp;
    VecGetSize( partialCurrent, &size );
    if ( size == numx )
    {
        //cout << " LegendrePoly: could be expanding a horizontal current " << endl;
        for ( integer i = 0; i <= maxOrder; i++ )
        {
            VecDot(Px[i], partialCurrent, &tmp);
            legCoefs[i] = tmp*wx[i];
        }
    }
    else if ( size == numy )
    {
        //cout << " LegendrePoly: could be expanding a vertical current " << endl;
        for ( integer i = 0; i <= maxOrder; i++ )
        {
            VecDot(Py[i], partialCurrent, &tmp);
            legCoefs[i] = tmp*wy[i];
        }
    }
    else
    {
        cout << " LegendrePolyERROR: input Vec matches no precomputed LP length! " << endl;
    }
    return;
}


//---------------------------------------------------------------------------//
/*!
 * \brief This is a simple factorial function.
 *
 */
scalar LegendrePoly::factorial( integer x ) 
{
     scalar fac = 1.0;
     for ( integer i=2; i<=x; i++ ) 
        fac *= i;
     return fac;
}


//---------------------------------------------------------------------------//
//                 end of LegendrePoly.cc
//---------------------------------------------------------------------------//

