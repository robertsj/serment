//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Diff2dProblem.hh
 * \author Jeremy Roberts
 * \date   10/26/2010
 * \brief  A class for handling input for diff2d problems.
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 171                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-11-19 20:36:52 -0500 (Sat, 19 Nov 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef DIFF2DPROBLEM_HH
#define DIFF2DPROBLEM_HH

#include <iostream>
#include <fstream>
#include <vector>
#include "petscvec.h"
#include "petscmat.h"
#include "petscsys.h"
#include "../linalg/typedefs.hh"
#include "Diff2dElement.hh"
#include "Diff2dInput.hh"

#include "utilities/SP.hh"



using namespace std;

//===========================================================================//
/*!
 * \class Diff2dProblem
 * \brief A class containing operators (i.e. matrices) for diff2d problems.
 *
 *  To be completed.
 *
 */
//===========================================================================//

class Diff2dProblem
{
  public:

        typedef util::SP<Diff2dProblem>   SP_problem;

        // multigroup diffusion operators (lhs operator)
        vector<Mat>             K;
        // inhomogeneous (fixed) source
        vector<Vec>             fixedSource;
        // homogeneous (i.e. fission) component = nuSigF
        vector<Vec>             nuSigF;  
        // the location and group dependent xi spectrum
        vector<Vec>             xiSpect;
        // the location and group dependent scattering cross-section
        vector<vector<Vec> >    sigScatter;  
        // the location and group dependent absorption cross-section
        vector<Vec>             sigAbs;
        // the location and group dependent scalar flux
        //vector<Vec>             phi;
        // x dimension (for easy manipulation)
        integer                 nx;
        // y dimension (for easy manipulation)
        integer                 ny;
        vector<scalar>          dx;
        vector<scalar>          dy;
        // the input
        Diff2dInput             *inp;
        // the element
        Diff2dElement           *el;
        scalar *cix;
        scalar *ciy;
        bool isRFSet;
        //--------------------------------------------------------------------

        Diff2dProblem(){};

        Diff2dProblem( Diff2dInput &input, integer elid );

        ~Diff2dProblem();
        // update RHS variables, only relevant for ptype 2
        void updateRHS( scalar k, int O, int G, int S);
        scalar getKeff(){ return keff; };
        void destroy();

    public:
        // response function specific quantities:
        // incident order
        int RO;
        // incident group
        int RG;
        // incident side
        int side;
        // keff for rhs
        scalar keff;
        //
        scalar *dv;

        // (TEMPORARY?) RESPONSE FUNCTION PLACEMENTS
        // everything in 1-D arrays for quick insertion into Petsc structures
        // in particular, R and RL represent "logically 2-D" arrays for use
        // in response matrix and leakage operators (as block components)
        scalar *R, *RL, *RF, *RA;
        scalar *Rold, *RLold, *RFold, *RAold;

        // k's associated with R's and Rold's
        scalar kNew, kOld;

        //--------------------------------------------------------------------
        //void setup();   // allocate petsc vecs and things
        void setK();    // build the matrix.  separate for multiple problems using one matrix set
        void setRHS();  // internal method for building the right hand side
        void setValues( scalar A[], integer idx[],  
                        integer k, integer lenk, integer g,
                        scalar dv, scalar cix_i, scalar ciy_j, 
                        integer m_C ); // why can't i pass these reference?
        // this creates the response function arrays to be used throughout
        void resizeR( integer m, scalar kk ) // i am temporary
        { 

            if (isRFSet==false) // first pass through
            {

                R  = new scalar[ (4*m) * (4*m) ]; Rold  = new scalar[ (4*m) * (4*m) ];
                RL = new scalar[ 4 * 4*m ];       RLold = new scalar[ 4 * 4*m ];
                RF = new scalar[ 4*m ];           RFold = new scalar[ 4*m ];
                RA = new scalar[ 4*m ];           RAold = new scalar[ 4*m ];
                isRFSet=true;
            }
            else
            {
                // switch-a-roo of the response function arrays
                kOld = kNew; kNew = kk; 
                scalar *tmpR;
                tmpR = R;  R  = Rold;  Rold  = tmpR;
                tmpR = RL; RL = RLold; RLold = tmpR;
                tmpR = RF; RF = RFold; RFold = tmpR;
                tmpR = RA; RA = RAold; RAold = tmpR;
            }
            clearR(m); // must reset given the way i fill them
           
        };
        // this resets the response function arrays to zero
        void clearR( integer m ) // i too am temporary
        { 
            for (integer i = 0; i < 4*m*4*m; i++)
                R[i]=0;
            for (integer i = 0; i < 4*4*m; i++)
                RL[i]=0;
            for (integer i = 0; i < 4*m; i++)
                RF[i]=0;
            for (integer i = 0; i < 4*m; i++)
                RA[i]=0;
        };
};

#endif // DIFF2DPROBLEM_HH

//---------------------------------------------------------------------------//
//                 end of Diff2dProblem.hh
//---------------------------------------------------------------------------//

