//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Connect2dCart.cc
 * \author Jeremy Roberts
 * \date   10/22/2010
 * \brief  Member definitions of concrete class Connect2dCart
 * \note   Copyright (C) 2010 Jeremy Roberts. \todo This is so ugly it makes me want to cry.
 */
//---------------------------------------------------------------------------//
// $Rev:: 171                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-11-19 20:36:52 -0500 (Sat, 19 Nov 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#include <iostream>
#include <math.h>
#include "petscvec.h"
#include "petscmat.h"
#include "LinAlg.hh"
#include "ConnectMatrix.hh"
#include "Connect2dCart.hh"

using namespace std;

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
Connect2dCart::Connect2dCart( GlobalInput::SP_globalinput in )
  : ConnectMatrix( &(*in) )
{
    cout << " CONSTRUCTING Connect2dCart " << endl;
    buildMe();
	return;
}

//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
Connect2dCart::~Connect2dCart()
{
    // nothing here right now, as M is destroyed explicitly in Release()
    return; 
}

//---------------------------------------------------------------------------//
/*!
 * \brief This builds the 2d Cartesian connectivity matrix.
 *
 * This routine takes uses the 2d array of element placements and produces
 * the associated connectivity matrix.
 * 
 * Suppose the problem is laid out as follows:
 *  \f[
 *  \left ( \begin{array}{ccc}
 *	  1 & 1 & 2  \\
 *	  1 & 1 & 2  \\
 *	  2 & 2 & 2
 *	\end{array} 
 *  \right ) \begin{array}{l}
       y \rightarrow \\
       x \\
       \downarrow  \end{array} \, ,
 *  \f]
 * where element "1" may be fuel and "2" a moderator.  Note, \f$ x \f$ and 
 * \f$ y \f$ take on
 * indices "i" and "j" so that the matrix looks transposed when written.  The
 * elements are ordered starting from the (physical) bottom left (i.e. the
 * upper left of the matrix above), and then along the \f$ x \f$ dimension for
 * one a single value of \f$ y \f$, i.e. physical row by physical row from
 * the bottom.  An extension to 3d would keep this structure for 
 * \f$ z \f$-planes, and then have an outer loop from bottom to top in the 
 * \f$ z \f$ direction. 
 *
 */
void Connect2dCart::buildMe()
{
    int I, J;  
    I = input->elements.size();
    J = input->elements[0].size();
    
    // size for the "zeroth order" connectivity matrix, mm
    integer mmSize = input->faces*input->numel;
    vector<vector<integer> > mm;
    mm.resize( mmSize );   
    for ( integer j = 0; j < mmSize; j++ )
        mm[j].resize( mmSize ); 

    for ( integer j = 0; j < mmSize; j++ )
        for ( integer i = 0; i < mmSize; i++ )
            mm[i][j] = 0;


    integer kk = 0; // size of little m, incremented for each reflective
    integer k = 0;
    integer nzrow = 0;
    for ( integer j = 0; j < J; j++ )        // "physical rows"
    {
        for ( integer i = 0; i < I; i++ )    // "physical columns"
        {
            if ( input->elements[i][j][0] >= 0 ) 
            {
                k++;  // I am the kth element; else, I'm just void so move on
          
		        // If there is somebody (>0) to my right, connect my face
		        // number 2 to their face number 1
		        if ( (i < I-1) and ( input->elements[i+1][j][0] >= 0 ) )
		        {
		            mm[ 4*(k-1)+1 ][ 4*k ] = 2;
		        }
		        // If there is somebody to my left, connect my face
		        // number 1 to their face number 2
		        if ( (i > 0) and ( input->elements[i-1][j][0] >= 0 ) )
		        {
		            mm[ 4*(k-1) ][ 4*(k-2)+1 ] = 3;
		        }
		        // If there is somebody to my top, connect my face 4 with
		        // their face 3.
		        if ( (j < J-1) and ( input->elements[i][j+1][0] >= 0 ) )
		        {
		            // nzrow counts elements left in my row and in the next
		            // row before my neighbor
		            nzrow = elSum( i, I, j, j+1 );
		            mm[ 4*(k-1)+3 ][ 4*(k-1+nzrow)+2 ] = 4;
		        }
		        // If somebody is to my bottom, my 3 to their 4
		        if ( (j > 0) and ( input->elements[i][j-1][0] >= 0 ) )
		        {
		            // nzrow counts elements left in my row and in the previous
		            // row before my neighbor
		            nzrow = elSum( i, I, j-1, j ); // had been j, j-1
		            mm[ 4*(k-1)+2 ][ 4*(k-1-nzrow)+3 ] = 5;
		        }

		        // boundaries
		        if ( (i == 0)   and ( input->bcl == 1) )      // left
		            mm[ 4*(k-1) ][ 4*(k-1) ]     = 1;
		        else
		            kk++;

		        if ( (i == I-1) and ( input->bcr == 1) )      // right
		            mm[ 4*(k-1)+1 ][ 4*(k-1)+1 ] = 1;
		        else
		            kk++;

		        if ( (j == 0)   and ( input->bcb == 1) )      // bot
		            mm[ 4*(k-1)+2 ][ 4*(k-1)+2 ] = 1;
		        else
		            kk++;

		        if ( (j == J-1) and ( input->bct == 1) )      // top
		            mm[ 4*(k-1)+3 ][ 4*(k-1)+3 ] = 1;     
		        else
		            kk++;   
            }
        }
    }

    // build little M
    littleMidx = new integer[kk]; 
    kk = 0;
    bool tmp;
    for ( integer j = 0; j < mmSize; j++ )      
    {
        tmp = false;
        for ( integer i = 0; i < mmSize; i++ )      
        {
            if ( mm[i][j] > 0 )   // if any value in column j is nonzero
                tmp = true;       // then there is either reflection or
        }                         // a connection, i.e. no leakage
        if ( tmp == false )       
        {                         // otherwise, tmp remains false 
            littleMidx[kk] = j;   // and the face corresponding to mm[j][j]
            kk++;                 // is adjacent to vacuum
        }
    }
    littleMsize = kk;

    // SETUP FOR REFLECTION:  Because the current responses are computed
    //   w/r to the INCOMING coordinates, the outgoing response is exactly
    //   the input of the neighbor cell; however, because reflection switches
    //   left-to-right to right-to-left, all the odd order Legendre moments
    //   must have their signs reversed.  Hence, the -1's here.
    // ***NOTE--not including angular explicitly here (since I don't actually
    //          know how it should be done)

    integer NumOrdGrps = ( input->spaceord+1)*( input->angleord+1)*input->numgroups;
    scalar refl[ NumOrdGrps ]; // a "logical" 2d array
    k = 0;
    for ( integer o = 0; o < input->spaceord+1; o++ )
    {
        for ( integer g = 0; g < input->numgroups; g++ )
        {
            refl[k] = pow( -1, o );
            k++;
        }
    }
    // For transmision, it's just a bunch of 1's
    scalar tran[ NumOrdGrps ];
    for ( integer i = 0; i < NumOrdGrps ; i++ )
        tran[i] = 1.0;
    // now, we need to insert into C.M a diagonal of "refl" wherever mm has 1 
    //  and a diagonal of "tran" wherever mm has 2.  The idea will be to 
    //  cycle through mm and build for each nonzero entry the index array
    //  for use with InsertVals .
    integer idxrow;
    integer idxcol;
    scalar val;
    //scalar time = MPI_Wtime();
    for ( integer j = 0; j < mmSize; j++ )
    {
        for ( integer i = 0; i < mmSize; i++ )
        {
            if ( mm[i][j] > 0 ) // build the index array
            {
                for (integer g = 0; g < NumOrdGrps; g++)
                {

                    idxrow = g + NumOrdGrps*i;
                    idxcol = g + NumOrdGrps*j;
                    if ( mm[i][j] == 1 )
                    {
                        val = refl[g];
                    }
                    else
                    {
                        val = tran[g];
                    }
                    insertVal( val, idxrow, idxcol );  
                }              
            }
            // else, move on
        }
    }
    // set M;
    checkReady();
}

//---------------------------------------------------------------------------//
/*!
 * \brief This counts non-zero elements between an element and a neighbor.
 *
 * Here, ii is the current row index, II is the number of rows, j1 is the
 * lower column, and j2 is the upper column.  Either j1 or j2 corresponds to 
 * to the current column index, depending if we are looking ahead for our 
 * neighbor to the top or looking below for our neighbor to the bottom.
 *
 */
integer Connect2dCart::elSum( integer ii, integer II, integer j1, integer j2 )
{
    integer sum = 0;
    for ( integer i = ii+1; i < II; i++ )
    {
        if ( input->elements[i][j1][0] >= 0 )
            sum++;
    }
    for ( integer i = 0; i <= ii; i++ )
    {
        if ( input->elements[i][j2][0] >= 0 )
            sum++;
    }
    return sum;
}

//---------------------------------------------------------------------------//
//                 end of Connect2dCart.cc
//---------------------------------------------------------------------------//

