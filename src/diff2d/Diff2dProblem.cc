//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Diff2dProblem.cc
 * \author Jeremy Roberts
 * \date   10/26/2010
 * \brief  Member definitions of class Diff2dProblem
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
#include "petscvec.h"
#include "petscmat.h"
#include "Diff2dElement.hh"
#include "Diff2dInput.hh"
#include "Diff2dProblem.hh"
#include "../linalg/typedefs.hh"

using namespace std;

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
Diff2dProblem::Diff2dProblem( Diff2dInput &input, integer elid )
{
    RO   = 0;
    RG   = 0;
    side = 0;
    keff = 1.0; kOld = -200.0; kNew = -100.0;
    isRFSet = false;

    //cout << " --- building the problem " << endl;

    inp = &input;
    el  = &inp->elements[elid];
    if ( inp->ptype == 2 )
        el->P.buildMe( el->nxfm, el->nyfm, inp->maxOrder );

    // matrix dimension
    integer n = el->nxfm * el->nyfm;
    K.resize( input.numg );

    for (integer g = 0; g < input.numg; ++g)
    {
        MatCreate( PETSC_COMM_SELF, &K[g] );
        MatSetType( K[g], MATAIJ );
        MatSetSizes( K[g], PETSC_DECIDE, PETSC_DECIDE, n, n ); 
        MatSetFromOptions( K[g] );
    }
    // build the matrix
    setK();
    //cout << " all done setting K " << endl; 

    return;
}

//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
Diff2dProblem::~Diff2dProblem()
{
    cout << " --- killing the problem " << endl;
    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief This function actually sets the K matrices
 *
 */
void Diff2dProblem::setK()
{

    dx.resize(el->nxfm);
    dy.resize(el->nyfm);
    dv  = new scalar[el->nxfm*el->nyfm];
    cix = new scalar[el->nxfm];
    ciy = new scalar[el->nyfm];

    // Build cix and ciy, which associate a coarse mesh index with each 
    // fine mesh index.  This makes array construction much cleaner.  Also
    // build the dx and dy vectors.
    integer ind_str, ind_end;
    integer j = 0;
    for( int i = 0; i < el->nxcm; ++i )
    {
        ind_str = j;
        ind_end = j + el->xfm[i];
        for( int gg = ind_str; gg < ind_end; ++gg )
        {
            dx[gg]  = ( el->xcm[i+1]-el->xcm[i] ) / el->xfm[i];
            cix[gg] = i;
        }
        j = 0;
        for( int gg = 0; gg <= i; ++gg )
            j = j + el->xfm[gg];
    }
    j = 0;
    for( integer i = 0; i < el->nycm; ++i )
    {
        ind_str = j;
        ind_end = j + el->yfm[i];
        for( integer gg = ind_str; gg < ind_end; ++gg )
        {
            dy[gg]  = ( el->ycm[i+1]-el->ycm[i] ) / el->yfm[i];
            ciy[gg] = i;
        }
        j = 0;
        for( integer gg = 0; gg <= i; ++gg )
            j = j + el->yfm[gg];
    }
    // the size of a matrix operator for a single group
    scalar lenk = el->nxfm * el->nyfm; 

    // allocate the various vectors of Vec's; this is all done here, because
    // any subsequent RHS-tweaking will happen after the matrix is done
    sigAbs.resize( inp->numg );
    nuSigF.resize( inp->numg );
    xiSpect.resize( inp->numg );
    sigScatter.resize( inp->numg ); 
    for ( integer i = 0; i < inp->numg; ++i )
        sigScatter[i].resize( inp->numg );
    fixedSource.resize( inp->numg );  

    // now allocate the consituent Vec's
    VecCreate( PETSC_COMM_SELF, &fixedSource[0] );
    VecSetSizes( fixedSource[0], PETSC_DECIDE, lenk );
    VecSetFromOptions( fixedSource[0] );
    VecDuplicate( fixedSource[0], &nuSigF[0] );
    VecDuplicate( fixedSource[0], &xiSpect[0] );
    VecDuplicate( fixedSource[0], &sigAbs[0] );  // only needed if for rf

    for ( integer i = 1; i < inp->numg; ++i )
    {
        VecDuplicate( fixedSource[0], &fixedSource[i] );
        VecDuplicate( fixedSource[0], &nuSigF[i]      );
        VecDuplicate( fixedSource[0], &xiSpect[i]     );
        VecDuplicate( fixedSource[0], &sigAbs[i]      );
    }

    for ( integer i = 0; i < inp->numg; ++i )
    {
        for ( integer j = 0; j < inp->numg; ++j )
        {
            VecDuplicate( fixedSource[0], &sigScatter[i][j] );
        }
    }
    
    // set Vec's of dx, dy, and dv, all part of Diff2dElement
    VecCreate( PETSC_COMM_WORLD, &el->dx );
    VecSetSizes( el->dx, PETSC_DECIDE, el->nxfm );
    VecSetFromOptions( el->dx );
    VecCreate( PETSC_COMM_WORLD, &el->dy );
    VecSetSizes( el->dy, PETSC_DECIDE, el->nxfm );
    VecSetFromOptions( el->dy );
    VecDuplicate( fixedSource[0], &el->dv  );
    integer k;
    for ( integer i = 0; i < el->nxfm; ++i )
    {
        VecSetValue( el->dx, i, dx[i], INSERT_VALUES );
        for ( integer j = 0; j < el->nyfm; ++j )
        {
                if ( i == 0 ) 
                {
                    VecSetValue( el->dy, j, dy[j], INSERT_VALUES );
                }
                k     = i + j*el->nxfm;
                dv[k] = dx[i] * dy[j];
                VecSetValue( el->dv, k, dv[k], INSERT_VALUES );
        }
    }
    VecAssemblyBegin( el->dx );
    VecAssemblyEnd( el->dx );
    VecAssemblyBegin( el->dy );
    VecAssemblyEnd( el->dy );
    VecAssemblyBegin( el->dv );
    VecAssemblyEnd( el->dv );


    //-----------------------------------------
    // BOUNDARY SOURCE: Here, I need to do two things.  One, if this is a stand
    // alone problem, then input has coarse-mesh boundary source definitions,
    // not implemented yet.  Here, those def's would be used to build the
    // needed alpha vectors directly.  If, however, this is a ptype=2 i.e.
    // response function problem, then the appropriate order LP in the el
    // object is placed in the appropriate group of the appropriate alphaXX.
    // 
    // Since the only thing that changes for a new side/RO/RG permutation is the
    // fixed source, it would be highly advantageous to have a separate source
    // construction routine.  Since this does will not change fundamentally 
    // the current structure, leave this as a TO DO.
    vector<vector<scalar> >  alphaL;
    vector<vector<scalar> >  alphaR;
    vector<vector<scalar> >  alphaB;
    vector<vector<scalar> >  alphaT;

    // alpha's ( boundary conditions )
    scalar one = 1.0;
    if ( el->bcl == 2 )
    {
        alphaL.resize( el->nyfm ); 
        for ( integer i = 0; i < el->nyfm; ++i )
        {
            alphaL[i].resize( inp->numg );
            fill( alphaL[i].begin(), alphaL[i].end(), one );
        }
    }
    if ( el->bcr == 2 )
    {
        alphaR.resize( el->nyfm ); 
        for ( integer i = 0; i < el->nyfm; ++i )
        {
            alphaR[i].resize( inp->numg );
            fill( alphaR[i].begin(), alphaR[i].end(), one );
        }
    }
    if ( el->bcb == 2 )
    {
        alphaB.resize( el->nxfm ); 
        for ( integer i = 0; i < el->nxfm; ++i )
        {
            alphaB[i].resize( inp->numg );
            fill( alphaB[i].begin(), alphaB[i].end(), one );
        }
    }
    if ( el->bct == 2 )
    {
        alphaT.resize( el->nxfm ); 
        for ( integer i = 0; i < el->nxfm; ++i )
        {
            alphaT[i].resize( inp->numg );
            fill( alphaT[i].begin(), alphaT[i].end(), one );
        }
    }

    // time to build!
    integer m_L, m_R, m_B, m_T, m_C;
    scalar  A[5];
    integer idx[4];
    scalar ntwo = -2.0;
    integer beta;

    //-----------------------------------------------------------------------
    // group loop over all coefficients

    for ( integer g = 0; g < inp->numg; ++g )
    {
        //-------------------------------------------------------------------
        // interior coefficients

        for( integer i = 1; i < el->nxfm-1; ++i )
        {
            for( integer j = 1; j < el->nyfm-1; ++j )
            {

                // materials
                m_L = el->mt[ cix[i-1] ][ ciy[j  ] ];
                m_R = el->mt[ cix[i+1] ][ ciy[j  ] ];
                m_B = el->mt[ cix[i  ] ][ ciy[j-1] ];
                m_T = el->mt[ cix[i  ] ][ ciy[j+1] ];
                m_C = el->mt[ cix[i  ] ][ ciy[j  ] ];

                // matrix row index
                k = i + j*el->nxfm;
                
                // define matrix values
               
                A[0] = ntwo*inp->dc[m_L][g]*inp->dc[m_C][g]*dy[j] / 
                      ( dx[i-1]*inp->dc[m_C][g] + dx[i]*inp->dc[m_L][g] );
              
                A[1] = ntwo*inp->dc[m_R][g]*inp->dc[m_C][g]*dy[j] / 
                      ( dx[i+1]*inp->dc[m_C][g] + dx[i]*inp->dc[m_R][g] );
                A[2] = ntwo*inp->dc[m_B][g]*inp->dc[m_C][g]*dx[i] / 
                      ( dy[j-1]*inp->dc[m_C][g] + dy[j]*inp->dc[m_B][g] );
                A[3] = ntwo*inp->dc[m_T][g]*inp->dc[m_C][g]*dx[i] /
                      ( dy[j+1]*inp->dc[m_C][g] + dy[j]*inp->dc[m_T][g] );
                A[4] = dv[k]*inp->sr[m_C][g] - ( A[0] + A[1] + A[2] + A[3] );

                // set values         
                idx[0] = 1; idx[1] = 1; idx[2] = 1; idx[3] = 1;
                setValues( A, idx, k, lenk, g, dv[k], cix[i], ciy[j], m_C );
            }
        }
     
        //-------------------------------------------------------------------
        // left and right edge coefficients

        for( integer j = 1; j < el->nyfm-1; ++j )
        {
            //---------------------------------------------------------------
            // left edge
            integer i = 0;
            k = i+j*el->nxfm;  

            if (el->bcl == 1)
                beta = 1;
            else
                beta = 0;
            // materials
            m_R = el->mt[ cix[i+1] ][ ciy[j  ] ];
            m_B = el->mt[ cix[i  ] ][ ciy[j-1] ];
            m_T = el->mt[ cix[i  ] ][ ciy[j+1] ];
            m_C = el->mt[ cix[i  ] ][ ciy[j  ] ];
            // coefficients
            A[1] = ntwo*inp->dc[m_R][g]*inp->dc[m_C][g]*dy[j] / 
                      ( dx[i+1]*inp->dc[m_C][g] + dx[i]*inp->dc[m_R][g] );
            A[2] = ntwo*inp->dc[m_B][g]*inp->dc[m_C][g]*dx[i] / 
                      ( dy[j-1]*inp->dc[m_C][g] + dy[j]*inp->dc[m_B][g] );
            A[3] = ntwo*inp->dc[m_T][g]*inp->dc[m_C][g]*dx[i] /
                      ( dy[j+1]*inp->dc[m_C][g] + dy[j]*inp->dc[m_T][g] );
            A[4] = dv[k]*inp->sr[m_C][g] - ( A[1] + A[2] + A[3] ) 
                  + inp->dc[m_C][g]*dy[j]*2.0*(1-beta) / 
                  ( 4.0*(1+beta)*inp->dc[m_C][g] + dx[i]*(1-beta) );
            // set values
            idx[0] = 0; idx[1] = 1; idx[2] = 1; idx[3] = 1;
            setValues( A, idx, k, lenk, g, dv[k], cix[i], ciy[j], m_C );

            //---------------------------------------------------------------
            // right edge
            i = el->nxfm - 1;
            k = i+j*el->nxfm;

            if (el->bcr == 1)
                beta = 1;
            else
                beta = 0;
            // materials
            m_L = el->mt[ cix[i-1] ][ ciy[j  ] ];
            m_B = el->mt[ cix[i  ] ][ ciy[j-1] ];
            m_T = el->mt[ cix[i  ] ][ ciy[j+1] ];
            m_C = el->mt[ cix[i  ] ][ ciy[j  ] ];
            // coefficients
            A[0] = ntwo*inp->dc[m_L][g]*inp->dc[m_C][g]*dy[j] / 
                      ( dx[i-1]*inp->dc[m_C][g] + dx[i]*inp->dc[m_L][g] );
            A[2] = ntwo*inp->dc[m_B][g]*inp->dc[m_C][g]*dx[i] / 
                      ( dy[j-1]*inp->dc[m_C][g] + dy[j]*inp->dc[m_B][g] );
            A[3] = ntwo*inp->dc[m_T][g]*inp->dc[m_C][g]*dx[i] /
                      ( dy[j+1]*inp->dc[m_C][g] + dy[j]*inp->dc[m_T][g] );
            A[4] = dv[k]*inp->sr[m_C][g] - ( A[0] + A[2] + A[3] ) 
                   + inp->dc[m_C][g]*dy[j]*2.0*(1-beta) / 
                   ( 4.0*(1+beta)*inp->dc[m_C][g] + dx[i]*(1-beta) );
            idx[0] = 1; idx[1] = 0; idx[2] = 1; idx[3] = 1;
            setValues( A, idx, k, lenk, g, dv[k], cix[i], ciy[j], m_C );

        } // end left and right coefficients

        //-------------------------------------------------------------------
        // top and bottom edge coefficients

        for( integer i = 1; i < el->nxfm-1; ++i )
        {
            //---------------------------------------------------------------
            // bottom edge
            integer j = 0;
            k = i+j*el->nxfm;

            if (el->bcb == 1)
                beta = 1;
            else
                beta = 0;
            // materials
            m_L = el->mt[ cix[i-1] ][ ciy[j  ] ];
            m_R = el->mt[ cix[i+1] ][ ciy[j  ] ];
            m_T = el->mt[ cix[i  ] ][ ciy[j+1] ];
            m_C = el->mt[ cix[i  ] ][ ciy[j  ] ];
            // coefficients
            A[0] = ntwo*inp->dc[m_L][g]*inp->dc[m_C][g]*dy[j] / 
                      ( dx[i-1]*inp->dc[m_C][g] + dx[i]*inp->dc[m_L][g] );
            A[1] = ntwo*inp->dc[m_R][g]*inp->dc[m_C][g]*dy[j] / 
                      ( dx[i+1]*inp->dc[m_C][g] + dx[i]*inp->dc[m_R][g] );
            A[3] = ntwo*inp->dc[m_T][g]*inp->dc[m_C][g]*dx[i] /
                      ( dy[j+1]*inp->dc[m_C][g] + dy[j]*inp->dc[m_T][g] );
            A[4] = dv[k]*inp->sr[m_C][g] - ( A[0] + A[1] + A[3] ) 
                  + inp->dc[m_C][g]*dx[i]*2.0*(1-beta) / 
                  ( 4.0*(1+beta)*inp->dc[m_C][g] + dy[j]*(1-beta) );
            // set values
            idx[0] = 1; idx[1] = 1; idx[2] = 0; idx[3] = 1;
            setValues( A, idx, k, lenk, g, dv[k], cix[i], ciy[j], m_C );

            //---------------------------------------------------------------
            // top edge
            j = el->nyfm - 1;
            k = i+j*el->nxfm; 

            if (el->bct == 1)
                beta = 1;
            else
                beta = 0;
            // materials
            m_L = el->mt[ cix[i-1] ][ ciy[j  ] ];
            m_R = el->mt[ cix[i+1] ][ ciy[j  ] ];
            m_B = el->mt[ cix[i  ] ][ ciy[j-1] ];
            m_C = el->mt[ cix[i  ] ][ ciy[j  ] ];
            // coefficients
            A[0] = ntwo*inp->dc[m_L][g]*inp->dc[m_C][g]*dy[j] / 
                      ( dx[i-1]*inp->dc[m_C][g] + dx[i]*inp->dc[m_L][g] );
            A[1] = ntwo*inp->dc[m_R][g]*inp->dc[m_C][g]*dy[j] / 
                      ( dx[i+1]*inp->dc[m_C][g] + dx[i]*inp->dc[m_R][g] );
            A[2] = ntwo*inp->dc[m_B][g]*inp->dc[m_C][g]*dx[i] / 
                      ( dy[j-1]*inp->dc[m_C][g] + dy[j]*inp->dc[m_B][g] );
            A[4] = dv[k]*inp->sr[m_C][g] - ( A[0] + A[1] + A[2] ) 
                  + inp->dc[m_C][g]*dx[i]*2.0*(1-beta) / 
                  ( 4.0*(1+beta)*inp->dc[m_C][g] + dy[j]*(1-beta) );
            idx[0] = 1; idx[1] = 1; idx[2] = 1; idx[3] = 0;
            setValues( A, idx, k, lenk, g, dv[k], cix[i], ciy[j], m_C );

        } // end top and bottom coefficients

        scalar beta1, beta2;
        //-------------------------------------------------------------------
        // bottom left corner

        integer i = 0;
        integer j = 0;
        k = i+j*el->nxfm; 

        if (el->bcl == 1)
            beta1 = 1;
        else
            beta1 = 0;
        if (el->bcb == 1)
            beta2 = 1;
        else
            beta2 = 0;
        // materials
        m_R = el->mt[ cix[i+1] ][ ciy[j  ] ];
        m_T = el->mt[ cix[i  ] ][ ciy[j+1] ];
        m_C = el->mt[ cix[i  ] ][ ciy[j  ] ];      
        // coefficients
        A[1] = ntwo*inp->dc[m_R][g]*inp->dc[m_C][g]*dy[j] / 
                  ( dx[i+1]*inp->dc[m_C][g] + dx[i]*inp->dc[m_R][g] );
        A[3] = ntwo*inp->dc[m_T][g]*inp->dc[m_C][g]*dx[i] /
                  ( dy[j+1]*inp->dc[m_C][g] + dy[j]*inp->dc[m_T][g] );
        A[4] = dv[k]*inp->sr[m_C][g] - (  A[1] + A[3] ) 
               + inp->dc[m_C][g]*dy[j]*2*(1-beta1) / 
                 ( 4*(1+beta1)*inp->dc[m_C][g] + dx[i]*(1-beta1) )
               + inp->dc[m_C][g]*dx[i]*2*(1-beta2) /
                 ( 4*(1+beta2)*inp->dc[m_C][g] + dy[j]*(1-beta2) );
        idx[0] = 0; idx[1] = 1; idx[2] = 0; idx[3] = 1;
        setValues( A, idx, k, lenk, g, dv[k], cix[i], ciy[j], m_C );

        //-------------------------------------------------------------------
        // top left corner

        i = 0;
        j = el->nyfm-1;
        k = i+j*el->nxfm;  

        if (el->bcl == 1)
            beta1 = 1;
        else
            beta1 = 0;
        if (el->bct == 1)
            beta2 = 1;
        else
            beta2 = 0;
        // materials
        m_R = el->mt[ cix[i+1] ][ ciy[j  ] ];
        m_B = el->mt[ cix[i  ] ][ ciy[j-1] ];
        m_C = el->mt[ cix[i  ] ][ ciy[j  ] ];      
        // coefficients
        A[1] = ntwo*inp->dc[m_R][g]*inp->dc[m_C][g]*dy[j] / 
                  ( dx[i+1]*inp->dc[m_C][g] + dx[i]*inp->dc[m_R][g] );
        A[2] = ntwo*inp->dc[m_B][g]*inp->dc[m_C][g]*dx[i] / 
                   ( dy[j-1]*inp->dc[m_C][g] + dy[j]*inp->dc[m_B][g] );
        A[4] = dv[k]*inp->sr[m_C][g] - (  A[1] + A[2] ) 
               + inp->dc[m_C][g]*dy[j]*2*(1-beta1) / 
                 ( 4*(1+beta1)*inp->dc[m_C][g] + dx[i]*(1-beta1) )
               + inp->dc[m_C][g]*dx[i]*2*(1-beta2) /
                 ( 4*(1+beta2)*inp->dc[m_C][g] + dy[j]*(1-beta2) );
        idx[0] = 0; idx[1] = 1; idx[2] = 1; idx[3] = 0;
        setValues( A, idx, k, lenk, g, dv[k], cix[i], ciy[j], m_C );

        //-------------------------------------------------------------------
        // bottom right corner

        i = el->nxfm-1; 
        j = 0;
        k = i+j*el->nxfm;

        if (el->bcr == 1)
            beta1 = 1;
        else
            beta1 = 0;
        if (el->bcb == 1)
            beta2 = 1;
        else
            beta2 = 0;
        // materials
        m_L = el->mt[ cix[i-1] ][ ciy[j  ] ];
        m_T = el->mt[ cix[i  ] ][ ciy[j+1] ];
        m_C = el->mt[ cix[i  ] ][ ciy[j  ] ];    
        // coefficients
        A[0] = ntwo*inp->dc[m_L][g]*inp->dc[m_C][g]*dy[j] / 
                  ( dx[i-1]*inp->dc[m_C][g] + dx[i]*inp->dc[m_L][g] );
        A[3] = ntwo*inp->dc[m_T][g]*inp->dc[m_C][g]*dx[i] /
                  ( dy[j+1]*inp->dc[m_C][g] + dy[j]*inp->dc[m_T][g] );
        A[4] = dv[k]*inp->sr[m_C][g] - (  A[0] + A[3] ) 
               + inp->dc[m_C][g]*dy[j]*2*(1-beta1) / 
                 ( 4*(1+beta1)*inp->dc[m_C][g] + dx[i]*(1-beta1) )
               + inp->dc[m_C][g]*dx[i]*2*(1-beta2) /
                 ( 4*(1+beta2)*inp->dc[m_C][g] + dy[j]*(1-beta2) );
        idx[0] = 1; idx[1] = 0; idx[2] = 0; idx[3] = 1;
        setValues( A, idx, k, lenk, g, dv[k], cix[i], ciy[j], m_C );

        //-------------------------------------------------------------------
        // top right corner

        i = el->nxfm-1;
        j = el->nyfm-1;
        k = i+j*el->nxfm;

        if (el->bcr == 1)
            beta1 = 1;
        else
            beta1 = 0;
        if (el->bct == 1)
            beta2 = 1;
        else
            beta2 = 0;
        // materials
        m_L = el->mt[ cix[i-1] ][ ciy[j  ] ];
        m_B = el->mt[ cix[i  ] ][ ciy[j-1] ];
        m_C = el->mt[ cix[i  ] ][ ciy[j  ] ];      
        // coefficients
        A[0] = ntwo*inp->dc[m_L][g]*inp->dc[m_C][g]*dy[j] / 
                  ( dx[i-1]*inp->dc[m_C][g] + dx[i]*inp->dc[m_L][g] );
        A[2] = ntwo*inp->dc[m_B][g]*inp->dc[m_C][g]*dx[i] / 
                   ( dy[j-1]*inp->dc[m_C][g] + dy[j]*inp->dc[m_B][g] );
        A[4] = dv[k]*inp->sr[m_C][g] - (  A[0] + A[2] ) 
               + inp->dc[m_C][g]*dy[j]*2*(1-beta1) / 
                 ( 4*(1+beta1)*inp->dc[m_C][g] + dx[i]*(1-beta1) )
               + inp->dc[m_C][g]*dx[i]*2*(1-beta2) /
                 ( 4*(1+beta2)*inp->dc[m_C][g] + dy[j]*(1-beta2) );
        idx[0] = 1; idx[1] = 0; idx[2] = 1; idx[3] = 0;
        setValues( A, idx, k, lenk, g, dv[k], cix[i], ciy[j], m_C );

        //-------------------------------------------------------------------
        // assemble

        MatAssemblyBegin( K[g], MAT_FINAL_ASSEMBLY );
        MatAssemblyEnd( K[g], MAT_FINAL_ASSEMBLY );

    }

    return;
}


//---------------------------------------------------------------------------//
/*!
 * \brief Set values of the matrix.
 *
 */
void Diff2dProblem::setValues( scalar A[], 
                               integer idx[], 
                               integer k, 
                               integer lenk, 
                               integer g,
                               scalar dv, 
                               scalar cix_i, 
                               scalar ciy_j, 
                               integer m_C )
{
    // define row index
    integer rownum = 1;
    integer rowidx[1];
    rowidx[0] = k;

    // define column indices
    integer colnum = 5;
    integer colidx[5];
    colidx[0] = rowidx[0]-1;        // L
    colidx[1] = rowidx[0]+1;        // R
    colidx[2] = rowidx[0]-el->nxfm;  // B
    colidx[3] = rowidx[0]+el->nxfm;  // T
    colidx[4] = rowidx[0];          // C
    if ( idx[0] == 0 )
        colidx[0] = -1;
    if ( idx[1] == 0 )
        colidx[1] = -1;
    if ( idx[2] == 0 )
        colidx[2] = -1;
    if ( idx[3] == 0 )
        colidx[3] = -1;      
    
    // const integer *row; const scalar *val; 
    integer row; scalar val;
    // set matrix values
    MatSetValues( K[g], rownum, rowidx, colnum, colidx, A, INSERT_VALUES );
    // scatter valuues
    row = k; 
    for ( integer gg = 0; gg < inp->numg; ++gg )
    {
        val = dv*inp->sc[m_C][g][gg]; 
        VecSetValue( sigScatter[g][gg], row, val, INSERT_VALUES );
    }
    // compute fixed source if needed
    if ( inp->ptype == 0 )
    {
        val = dv*el->src[cix_i][ciy_j][g]; 
        VecSetValue( fixedSource[g], row, val, INSERT_VALUES );
    }
    val = dv*inp->ns[m_C][g];
    VecSetValue( nuSigF[g], row, val, INSERT_VALUES );
    val = inp->xi[m_C][g];
    VecSetValue( xiSpect[g], row, val, INSERT_VALUES );
    val = dv*inp->ab[m_C][g];
    VecSetValue( sigAbs[g], row, val, INSERT_VALUES );
    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief This function sets things on the RHS affected by updates
 *
 * The purpose of a separate fixedSource evaluator is to separate source 
 * terms updated during rf-generation from K and other quantities that are 
 * static for a given geometry/material definition.  This only includes the
 * fixedSource---all other quantities are static. By separating, we can
 * avoid the cost of re-constructing things.  For now, I'm setting it up only
 * for response function generation, i.e. I'm not accounting for stand alone
 * problems with a boundary source.
 *
 */
void Diff2dProblem::setRHS()
{
    scalar  *alphaV; // alphaV and alphaH point to the array of whichever
    scalar  *alphaH; //   Legendre polynomial is to be used as a source
    integer *indx;   // index for placement into petsc Vec
    scalar  *bsrc;   // actual values for the right hand side
    integer num;     // number of non zero values to be inserted
    scalar  zero = 0.0;

    // alpha's ( boundary conditions ), a scalar array for insertion
    // into a petsc Vec.  Get alpha from LegendrePoly.
   // cout << " NEW RHS, RO=" << RO << " RG="<< RG  << " side=" << side << endl;

    if ( side < 2 )  // LEFT OR RIGHT
        VecGetArray( el->P.Py[ RO ], &alphaV );
    else             // BOTTOM OR TOP
        VecGetArray( el->P.Px[ RO ], &alphaH );

    if ( side < 2 )
        num = el->nyfm;
    else
        num = el->nxfm;

    indx = new integer[ num ];
    bsrc = new scalar[ num ];

    integer m_C;
    integer idx;

    //-----------------------------------------------------------------------
    // group loop over all coefficients
    for ( integer g = 0; g < inp->numg; ++g )
    {

        VecSet( fixedSource[g], zero );
        VecAssemblyBegin( fixedSource[g] );
        VecAssemblyEnd( fixedSource[g] );
        if ( g == RG )
        {
            //----------------------------------------------------------------
            // left and right edge coefficients
            idx = 0;
            for( integer j = 1; j < el->nyfm-1; ++j )
            {
                //------------------------------------------------------------
                // left edge
                integer i = 0;
                m_C = el->mt[ cix[i  ] ][ ciy[j  ] ];
                if ( side == 0 )
                {
                    bsrc[idx] = 8.*inp->dc[m_C][g]*dy[j]*alphaV[el->nyfm-1-j] /
                                (4.*inp->dc[m_C][g]+dx[i]);
                    indx[idx] = i+j*el->nxfm;
                    idx++;
                }
                //------------------------------------------------------------
                // right edge
                i = el->nxfm - 1;
                m_C = el->mt[ cix[i  ] ][ ciy[j  ] ];
                if ( side == 1 )
                {
                    bsrc[idx] = 8.*inp->dc[m_C][g]*dy[j]*alphaV[j] / 
                               (4.*inp->dc[m_C][g]+dx[i]);
                    indx[idx] = i+j*el->nxfm;
                    idx++;
                }
            } 

            //----------------------------------------------------------------
            // top and bottom edge coefficients
            for( integer i = 1; i < el->nxfm-1; ++i )
            {
                //------------------------------------------------------------
                // bottom edge
                integer j = 0;
                m_C = el->mt[ cix[i  ] ][ ciy[j  ] ];
                if ( side == 2 )
                {
                    bsrc[idx] = 8.*inp->dc[m_C][g]*dx[i]*alphaH[i] / 
                                (4.*inp->dc[m_C][g]+dy[j]);
                    indx[idx] = i+j*el->nxfm;
                    idx++;
                }
                //------------------------------------------------------------
                // top edge
                j = el->nyfm - 1;
                m_C = el->mt[ cix[i  ] ][ ciy[j  ] ];
                if ( side == 3 )
                {
                    bsrc[idx] = 8.*inp->dc[m_C][g]*dx[i]*alphaH[el->nxfm-1-i] / 
                                (4.*inp->dc[m_C][g]+dy[j]);
                    indx[idx] = i+j*el->nxfm;
                    idx++;
                }
            } 

            //----------------------------------------------------------------
            // bottom left corner
            integer i = 0;
            integer j = 0;
            m_C = el->mt[ cix[i  ] ][ ciy[j  ] ];      
            if ( side == 0 )
            {
                bsrc[idx] = 8.*inp->dc[m_C][g]*dy[j]*alphaV[el->nyfm-1-j] / 
                            (4.*inp->dc[m_C][g]+dx[i]);
                indx[idx] = i+j*el->nxfm;
                idx++;
            }
            if ( side == 2 )
            {
                bsrc[idx] = 8.*inp->dc[m_C][g]*dx[i]*alphaH[i] / 
                            (4.*inp->dc[m_C][g]+dy[j]);
                indx[idx] = i+j*el->nxfm;
                idx++;
            }

            //----------------------------------------------------------------
            // top left corner
            i = 0;
            j = el->nyfm-1;
            m_C = el->mt[ cix[i  ] ][ ciy[j  ] ];      
            if ( side == 0 )
            {
                bsrc[idx] = 8.*inp->dc[m_C][g]*dy[j]*alphaV[el->nyfm-1-j] / 
                            (4.*inp->dc[m_C][g]+dx[i]);
                indx[idx] = i+j*el->nxfm;
                idx++;
            }
            if ( side == 3 )
            {
                bsrc[idx] = 8.*inp->dc[m_C][g]*dx[i]*alphaH[el->nxfm-1-i] / 
                            (4.*inp->dc[m_C][g]+dy[j]);
                indx[idx] = i+j*el->nxfm;
                idx++;
            }

            //----------------------------------------------------------------
            // bottom right corner
            i = el->nxfm-1; 
            j = 0;
            m_C = el->mt[ cix[i  ] ][ ciy[j  ] ];    
            if ( side == 1 )
            {
                bsrc[idx] = 8.*inp->dc[m_C][g]*dy[j]*alphaV[j] / 
                            (4.*inp->dc[m_C][g]+dx[i]);
                indx[idx] = i+j*el->nxfm;
                idx++;
            }
            if ( side == 2 )
            {
                bsrc[idx] = 8.*inp->dc[m_C][g]*dx[i]*alphaH[i] / 
                            (4.*inp->dc[m_C][g]+dy[j]);
                indx[idx] = i+j*el->nxfm;
                idx++;
            }

            //----------------------------------------------------------------
            // top right corner
            i = el->nxfm-1;
            j = el->nyfm-1;
            m_C = el->mt[ cix[i  ] ][ ciy[j  ] ];      
            if ( side == 1 )
            {
                bsrc[idx] = 8.*inp->dc[m_C][g]*dy[j]*alphaV[j] / 
                            (4.*inp->dc[m_C][g]+dx[i]);
                indx[idx] = i+j*el->nxfm;
                idx++;
            }
            if ( side == 3 )
            {
                bsrc[idx] = 8.*inp->dc[m_C][g]*dx[i]*alphaH[ el->nxfm-1-i] / 
                            (4.*inp->dc[m_C][g]+dy[j]);
                indx[idx] = i+j*el->nxfm;
                idx++;
            }
            
            //----------------------------------------------------------------
            // insert values
            VecSetValues( fixedSource[g], num, indx, bsrc, INSERT_VALUES );
            VecAssemblyBegin( fixedSource[g] );
            VecAssemblyEnd( fixedSource[g] );
    
        } // end if
    } // end for over groups

    // restore array
    if ( side < 2 ) // LEFT OR RIGHT
        VecRestoreArray( el->P.Py[ RO ], &alphaV );
    else            // BOTTOM OR TOP
        VecRestoreArray( el->P.Px[ RO ], &alphaH );

    // should I delete bsrc and indx??
    delete [] bsrc;
    delete [] indx;
    bsrc = NULL;
    indx = NULL;

    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Update the RHS with new order, group, and side
 *
 */
void Diff2dProblem::updateRHS( scalar k, int O, int G, int S )
{
    keff = k;
    RO = O;
    RG = G;
    side = S;
    setRHS();
    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destroy the K matrices ... finish me!
 *
 */
void Diff2dProblem::destroy()
{
    cout << " problem.destroy() ! " << endl;
    integer imax = K.size();
    integer jmax = sigScatter[0].size();
    for (integer i = 0; i < imax; ++i)
    {
        MatDestroy( &K[i] );
        VecDestroy( &fixedSource[i] );
        VecDestroy( &nuSigF[i] );
        VecDestroy( &xiSpect[i] );
        VecDestroy( &sigAbs[i] );
        //VecDestroy( phi[i] );
        for (integer j = 0; j < jmax; ++j)
        {
            VecDestroy( &sigScatter[i][j] );
        }
    }
    delete [] dv;
    delete [] cix;
    delete [] ciy;
    if ( isRFSet )
    {
        delete [] R;
        delete [] RL;
        delete [] RF;
        delete [] RA;
    }
    el->destroy();
    return;
}

//---------------------------------------------------------------------------//
//                 end of Diff2dProblem.cc
//---------------------------------------------------------------------------//

