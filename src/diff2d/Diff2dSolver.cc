//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Diff2dSolver.cc
 * \author Jeremy Roberts
 * \date   10/26/2010
 * \brief  Member definitions of class Diff2dSolver
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
#include "petscmat.h"
#include "petscksp.h"
#include "petsc.h"
#include "Diff2dInput.hh"
#include "Diff2dProblem.hh"
#include "Diff2dSolver.hh"
#include "../linalg/typedefs.hh"

using namespace std;

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
Diff2dSolver::Diff2dSolver(Diff2dInput &input, Diff2dProblem &problem, int el)
{
    //cout << " --- solving the problem " << endl;
    in = &input;   // pointer to input
    pr = &problem; // pointer to problem
    elid = el;
    nrow = in->elements[el].nxfm * in->elements[el].nyfm;
    ptype = in->ptype;
    // get the flux vectors ready
    phi.resize( in->numg );
    phi0.resize( in->numg );
    VecDuplicate( pr->nuSigF[0], &scsrc );
    VecDuplicate( pr->nuSigF[0], &tmp   );
    VecDuplicate( pr->nuSigF[0], &fisrc );
    VecDuplicate( pr->nuSigF[0], &fisrc0 );
    scalar zero = 0.0;
    scalar one = 1.0;
    for ( integer g = 0; g < in->numg; ++g )
    {
        VecDuplicate( pr->nuSigF[0], &phi[g] );
        VecDuplicate( pr->nuSigF[0], &phi0[g] );    
        VecSet( phi[g], zero );
        VecSet( phi0[g], zero );
        //VecAssemblyBegin( phi[g] ); VecAssemblyEnd( phi[g] );
    }
    VecSet( fisrc, one );
    VecAssemblyBegin( tmp ); VecAssemblyEnd( tmp );
    
    // initialize the solver context
    KSPCreate( PETSC_COMM_SELF, &ksp );
    // I'm using defaults here; they should be set by user input
    KSPSetTolerances(ksp,1e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT );
    KSPSetType(ksp,KSPGMRES);
    KSPGMRESSetRestart(ksp,10);

    //KSPGetPC(  ksp, &prec ); 
    //PCSetType( prec, PCSOR );
    //PCILUSetLevels( prec, 1 );
    //PCILUSetUseDropTolerance(prec,0.0001,0.05,5);

  //  KSPSetFromOptions( ksp );

    return;
}

//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
Diff2dSolver::~Diff2dSolver()
{
    cout << " --- killing the solver " << endl;

    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief This function drives the solution for a single problem.
 *
 */
void Diff2dSolver::solve()
{

    //pr->setK(); // set the matrices before solving.
    if ( ptype != 1 )
    {
        this->solveFixedMult();
    }
    else if ( ptype == 1 )
    {
        this->solveEigenvalue();
    }

    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief This function solves a fixed source problem with multiplication.
 *
 * The standard source iteration technique is used on the groups.  Upscattering
 * is accounted for (except into the first group, but this is a quick change;
 * it is omitted to make one group problems easier). 
 *
 */
void Diff2dSolver::solveFixedMult()
{
    //cout << " solving a fixed source problem.... " << in->numg << endl;

    scalar zero = 0.0;   // I probably use these all over...where is a
    scalar none = -1.0;  // good place for the definition?
    scalar one  = 1.0;

    VecSet( fisrc, zero ); // start with a *zero* guess and let build up
                           // probably better ways...think about it
    scalar k = pr->getKeff();  

    scalar norm = 1.0;
    integer it = 0;
    scalar errs;
    //cout << "----------------------------------------------" << endl;

    scalar tsolves = 0.0, tfiss=0.0, tcopy=0.0, ttmp=0.0, tsol2 = 0.0, ttmp2 = 0.0;

    do
    {

        // add an iteration
        it++;
       
        // begin group loop
        for ( integer g = 0; g < in->numg; ++g )
        {

            // set group operator
            KSPSetOperators( ksp, pr->K[g], pr->K[g], SAME_NONZERO_PATTERN );
            VecSet( tmp, zero );

            // special group 1 treatment (no upscatter...should change)
            if ( g == 0 )
            {
                // tmp = chi * fission source
                VecPointwiseMult( tmp, pr->xiSpect[g], fisrc );
                // tmp =  one/keff*tmp + fixed source
                VecAYPX( tmp, one/k, pr->fixedSource[g] );
                KSPSolve( ksp, tmp, phi[g] );
            }
            else
            {
                VecSet( scsrc, zero );
                for ( integer gg = 0; gg < in->numg; ++gg )
                {
                    VecPointwiseMult( tmp, phi[gg], pr->sigScatter[g][gg] );
                    VecAXPY( scsrc, one, tmp );
                }
                // tmp = chi * fission source
                VecPointwiseMult( tmp, pr->xiSpect[g], fisrc );
                // tmp = tmp/keff + fixed source
                VecAYPX( tmp, one/k, pr->fixedSource[g] );
                // tmp = tmp + scatter source = fission/k + fixed + scatter
                VecAYPX( tmp, one, scsrc );
                KSPSolve( ksp, tmp, phi[g] );
            }
        } // end group loop
    //    tsolves = tsolves + ( MPI_Wtime()-ttmp );

    //    ttmp = MPI_Wtime();
        // reset fission source
        VecSet( fisrc, zero );
        for ( integer g = 0; g < in->numg; ++g )
        {
            VecPointwiseMult( tmp, pr->nuSigF[g], phi[g] );
            VecAYPX( fisrc, one, tmp );
        }

        // find max phi error
        norm = 0.0;
        for ( integer g = 0; g < in->numg; ++g )
        {
            VecWAXPY( tmp, none, phi[g], phi0[g] );
            VecNorm( tmp, NORM_1, &errs );
            if ( errs > norm )
                norm = errs;
        }

        // save flux
        for ( integer g = 0; g < in->numg; ++g )
            VecCopy( phi[g], phi0[g] );
        
    }
    while ( ( norm > in->epss ) and ( it < in->maxit ) );

    return;
}


//---------------------------------------------------------------------------//
/*!
 * \brief This routine solves an eigenvalue problem.
 *
 * It uses a standard inner iteration over the groups and an outer iteration 
 * to update $k_{e\!f\!f}$.  In the future, possible acceleration techniques
 * could be used, including Chebyschev acceleration or Krylov subspace 
 * method such as the implicitly-restarted Arnoldi algorithm.
 *
 */
void Diff2dSolver::solveEigenvalue()
{
    cout << " solving an eigenvalue problem.... " << in->numg << endl;

    scalar zero = 0.0;   // I probably use these all over...where is a
    scalar none = -1.0;  // good place for the definition?
    scalar one  = 1.0;
    VecSet( fisrc, one );
    keff = 1.0;          // an initial guess should be given by user!!
    scalar norm = 1.0;
    integer it = 0;
    scalar keff0, errk, errs, gain, loss;
    cout << "----------------------------------------------" << endl;
    do
    {

        // add an iteration
        it++;

        // begin group loop
        for ( integer g = 0; g < in->numg; ++g )
        {

            // set group operator
            KSPSetOperators( ksp, pr->K[g], pr->K[g], DIFFERENT_NONZERO_PATTERN );
            VecSet( tmp, zero ); //VecAssemblyBegin( tmp ); VecAssemblyEnd( tmp );

            // special group 1 treatment (no scatter...should change)
            if ( g == 0 )
            {
                VecPointwiseMult( tmp, pr->xiSpect[g], fisrc );
                VecScale( tmp, one/keff );
                KSPSolve( ksp, tmp, phi[g] );
            }
            else
            {
                VecSet( scsrc, zero );
                for ( integer gg = 0; gg < in->numg; ++gg )
                {
                    VecPointwiseMult( tmp, phi[gg], pr->sigScatter[g][gg] );
                    VecAXPY( scsrc, one, tmp );
                }
                VecPointwiseMult( tmp, pr->xiSpect[g], fisrc );
                VecScale( tmp, one/keff );
                VecAYPX( tmp, one, scsrc );
                KSPSolve( ksp, tmp, phi[g] );
            }
        } // end group loop

        // save fission source and keff 
        VecCopy( fisrc, fisrc0 );
        keff0 = keff;
        // reset fission source
        VecSet( fisrc, zero );
        for ( integer g = 0; g < in->numg; ++g )
        {
            VecPointwiseMult( tmp, pr->nuSigF[g], phi[g] );
            VecAYPX( fisrc, one, tmp );
        }
        VecDot( in->elements[elid].dv, fisrc, &gain );
        VecDot( in->elements[elid].dv, fisrc0, &loss );

        keff = keff0*gain/loss;
        // find max fisrc error
        VecWAXPY( tmp, none, fisrc, fisrc0 );
        VecNorm( tmp, NORM_1, &errs );
        errk = abs( (keff-keff0)/keff );
        
        if ( in->printout == 1 and it % 10 == 0 )
        {
            printf (" iteration: %5i errs:  %8.6g  errk %8.6f keff %12.9f \n", 
                it, errs, errk, keff );

        }
    }
    while ( ( norm > in->epss ) and ( it < in->maxit ) );

    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief This routine computes the response functions.
 *
 * The current responses are computed.  The absorption, leakage, and fission
 * responses are computed.  Add some math here, maybe.  Remember, odd orders
 * are reversed for the incident side and for its right hand neighbor.
 *
 * For the current responses, note that the incident basis functions are 
 * always oriented with respect to the incident face.  That means for the 
 * top and left edges, the basis functions are reversed when the right hand
 * side is constructed.  That is accounted for here where the alpha terms are
 * indexed in reverse for those edges.  Similarly, the right and bottom 
 * current responses are constructed in reverse.  This is because the responses
 * must be expanded with respect to their face "looking outward" into the
 * next element.  Since top and left are already reversed, we expand as normal,
 * but for right and bottom, we must explicitly account for that.
 *
 */
void Diff2dSolver::compRespFct()
{
    
    integer k;
    integer m_C;
    scalar *phiArr;
    scalar phiPrimeOld;

    // Here, we pull out the Legendre Polynomial we used as an incident
    // current.  Below, we use it element-wise with a bit of logic. 
    scalar  *alphaV; // vertical LP
    scalar  *alphaH; // horz LP
    if ( pr->side == 0 or pr->side == 1 )
        VecGetArray( pr->el->P.Py[ pr->RO ], &alphaV );
    else
        VecGetArray( pr->el->P.Px[ pr->RO ], &alphaH );

    // These are Petsc Vec's into which we will place directly the CR arrays
    // without using new memory.  Doing so allows us to use Petsc's dot product
    // through use of the expandCur method of LegendrePoly.
    Vec crV, crH;
    VecDuplicate( pr->el->P.Px[ pr->RO ], &crH );
    VecDuplicate( pr->el->P.Py[ pr->RO ], &crV );

    // Temporary exiting partial currents
    scalar CRL[ pr->el->nyfm ];
    scalar CRR[ pr->el->nyfm ];
    scalar CRB[ pr->el->nxfm ];
    scalar CRT[ pr->el->nxfm ];
    scalar LL=0, LR=0, LB=0, LT=0;
    // Temporary array to hold expansion coefficients
    scalar FF[ in->maxOrder+1 ];

    // helpful quantities for indexing
    integer space  = (in->maxOrder+1)*in->numg;
    integer column =  pr->RO*in->numg+pr->RG + pr->side*space;
    integer ii, jj;

    // temporary for getting dot production for RF and RA
    scalar fisabs; 

    for ( integer g = 0; g < in->numg; g++ )
    {
        // get the array from Vec
        VecGetArray( phi[ g ], &phiArr );
        for ( integer j = 0; j < pr->el->nyfm; j++ )
        {
            //----------------------------------------------------------------
            // compute the left responses
            integer i = 0;
            k         = i+j*pr->el->nxfm;  
            CRL[j]    = 0.0;
            m_C       = pr->el->mt[ pr->cix[i  ] ][ pr->ciy[j  ] ]; // prolly want function in el or prob to give mt(i,j) for fine meshes
            phiPrimeOld = 0.0;
            if ( pr->side == 0 and g == pr->RG )
                phiPrimeOld = -8.0*alphaV[pr->el->nyfm-1-j];
            phiPrimeOld = ( 2.0 * phiArr[k] + phiPrimeOld ) /  
                                        ( 4.0*in->dc[m_C][g] + pr->dx[i] );
            // leakage
            LL          = LL + ( in->dc[m_C][g] * phiPrimeOld ) * pr->dy[j];     
            // current
            CRL[j]      = 0.25*phiArr[k] + 
                         ( 0.5*in->dc[m_C][g] - 0.125*pr->dx[i] ) * phiPrimeOld;
            //----------------------------------------------------------------
            // compute the right responses
            jj          = pr->el->nyfm-1-j; // this reverses placement
            i           = pr->el->nxfm-1;
            k           = i+j*pr->el->nxfm; 
            CRR[jj]     = 0.0; 
            m_C         = pr->el->mt[ pr->cix[i  ] ][ pr->ciy[j  ] ]; 
            phiPrimeOld = 0.0;
            if ( pr->side == 1 and g == pr->RG )
                phiPrimeOld = 8.0*alphaV[j];
            phiPrimeOld = (-2.0 * phiArr[k] + phiPrimeOld ) /  
                                       ( 4.0*in->dc[m_C][g] + pr->dx[i] );
            // leakage
            LR          = LR - ( in->dc[m_C][g] * phiPrimeOld ) * pr->dy[j];     
            // current
            CRR[jj]     = 0.25*phiArr[k] - 
                         ( 0.5*in->dc[m_C][g] - 0.125*pr->dx[i] ) * phiPrimeOld;
        }
        for ( integer i = 0; i < pr->el->nxfm; i++ )
        {
            //----------------------------------------------------------------
            // compute the bottom responses
            ii          = pr->el->nxfm-1-i; // this reverses
            integer j   = 0;
            k           = i+j*pr->el->nxfm;  
            m_C         = pr->el->mt[ pr->cix[i  ] ][ pr->ciy[j  ] ]; 
            CRB[ii]     = 0.0; 
            phiPrimeOld = 0.0;
            if ( pr->side == 2 and g == pr->RG )
                phiPrimeOld = -8.0*alphaH[i];
            phiPrimeOld = ( 2.0 * phiArr[k] + phiPrimeOld ) / 
                                        ( 4.0*in->dc[m_C][g] + pr->dy[j] );
            // leakage
            LB          = LB + ( in->dc[m_C][g] * phiPrimeOld ) * pr->dx[i];
            // current
            CRB[ii]     = 0.25*phiArr[k] + 
                         ( 0.5*in->dc[m_C][g] - 0.125*pr->dy[j] ) * phiPrimeOld;
            //----------------------------------------------------------------
            // compute the top response
            j           = pr->el->nyfm-1;
            k           = i+j*pr->el->nxfm;  
            
            m_C         = pr->el->mt[ pr->cix[i  ] ][ pr->ciy[j  ] ]; 
            CRT[i]      = 0.0;
            phiPrimeOld = 0.0;
            if ( pr->side == 3 and g == pr->RG )
                phiPrimeOld = 8.0*alphaH[pr->el->nxfm-1-i];
            phiPrimeOld = (-2.0 * phiArr[k] + phiPrimeOld ) /  
                                       ( 4.0*in->dc[m_C][g] + pr->dy[j] );
            // leakage
            LT          = LT - ( in->dc[m_C][g] * phiPrimeOld ) * pr->dx[i];     
            // current
            CRT[i]      = 0.25*phiArr[k] - 
                         ( 0.5*in->dc[m_C][g] - 0.125*pr->dy[j] ) * phiPrimeOld;
        }

        // put CRL in crV and expand;
        VecPlaceArray( crV, CRL );
        pr->el->P.expandCur( crV, FF );
        VecResetArray( crV );
        for ( integer o = 0; o < in->maxOrder+1; o++ )
            pr->R[ (g + o*in->numg + 0*space) + column*4*space ] = FF[o]; 
        //  pr->R[ g + o*in->numg + 0*space ][ column ] = FF[o]; 

        // put CRR in crV and expand;
        VecPlaceArray( crV, CRR );
        pr->el->P.expandCur( crV, FF );
        VecResetArray( crV );
        for ( integer o = 0; o < in->maxOrder+1; o++ )
            pr->R[ (g + o*in->numg + 1*space) + column*4*space ] = FF[o]; 
        // pr->R[ g + o*in->numg + 1*space ][ column ] = FF[o]; 

        // put CRB in crV and expand;
        VecPlaceArray( crH, CRB );
        pr->el->P.expandCur( crH, FF );
        VecResetArray( crH );
        for ( integer o = 0; o < in->maxOrder+1; o++ )
            pr->R[ (g + o*in->numg + 2*space) + column*4*space ] = FF[o]; 
        //  pr->R[ g + o*in->numg + 2*space ][ column ] = FF[o]; 

        // put CRT in crV and expand;
        VecPlaceArray( crH, CRT );
        pr->el->P.expandCur( crH, FF );
        VecResetArray( crH );
        for ( integer o = 0; o < in->maxOrder+1; o++ )
            pr->R[ (g + o*in->numg + 3*space) + column*4*space ] = FF[o]; 
        //  pr->R[ g + o*in->numg + 3*space ][ column ] = FF[o]; 

        // put the array back so we can dot product fission and absorption
        VecRestoreArray( phi[ g ], &phiArr );

        // FISSION
        VecDot( phi[ g ], pr->nuSigF[ g ], &fisabs );
        pr->RF[ column ] = fisabs + pr->RF[ column ];
        // ABSORPTION
        VecDot( phi[ g ], pr->sigAbs[ g ], &fisabs );
        pr->RA[ column ] = fisabs + pr->RA[ column ];        

    } // end group loop

    // place leakage responses
    pr->RL[0 + 4*column] = LL;
    pr->RL[1 + 4*column] = LR;
    pr->RL[2 + 4*column] = LB;
    pr->RL[3 + 4*column] = LT;
//    pr->RL[0][column] = LL;
//    pr->RL[1][column] = LR;
//    pr->RL[2][column] = LB;
//    pr->RL[3][column] = LT;

    if ( pr->side == 0 or pr->side == 1 )
        VecRestoreArray( pr->el->P.Py[ pr->RO ], &alphaV );
    else
        VecRestoreArray( pr->el->P.Px[ pr->RO ], &alphaH );

}

//---------------------------------------------------------------------------//
/*!
 * \brief Destroy the K matrices ... finish me!
 *
 */
void Diff2dSolver::destroy()
{
    cout << " solver.destroy() ! " << endl;
   
    for (integer i = 0; i < in->numg; ++i)
    {
        VecDestroy( &phi[i] );
        VecDestroy( &phi0[i] );
    }
    VecDestroy( &scsrc );
    VecDestroy( &fisrc );
    VecDestroy( &tmp );
    VecDestroy( &fisrc0 );

    return;
}

//---------------------------------------------------------------------------//
//                 end of Diff2dSolver.cc
//---------------------------------------------------------------------------//

