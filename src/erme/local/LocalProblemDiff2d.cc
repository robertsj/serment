//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LocalProblemDiff2d.cc
 * \author Jeremy Roberts
 * \date   11/24/2010
 * \brief  Member definitions of class LocalProblemDiff2d
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 140                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-09-14 12:53:40 -0400 (Wed, 14 Sep 2011) $:Date of last commit
//---------------------------------------------------------------------------//
#include <iostream>
#include <iomanip>
#include "LocalProblem.hh"
#include "LocalProblemDiff2d.hh"
#include "ResponseFunction.hh"
#include "ResponseFunctionDiffusion.hh"
#include "linalg/typedefs.hh"
#include "diff2d/Diff2dInput.hh"
#include "diff2d/Diff2dProblem.hh"
#include "diff2d/Diff2dSolver.hh"
#include "diff2d/Diff2dOutput.hh"
#include "petscmat.h"

using namespace std;
//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
LocalProblemDiff2d::LocalProblemDiff2d( string in )
    : LocalProblem(in)
{

    // process local input
    input.readInput( (char *)in.c_str() );

    keffCurrent = -100.0;
    keffOld     = -200.0;

    // create the local problem and solver vectors
    problem.resize( input.numel );
    solver.resize( input.numel );

    // initialize all the problems
    for (integer i = 0; i < input.numel; ++i)
    {
        problem[i] = new Diff2dProblem( input, i );
        solver[i]  = new Diff2dSolver( input, *problem[i], i ); 
    }

    // initialize the RFs
    RFs = new ResponseFunction*[input.numel];

}

//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
LocalProblemDiff2d::~LocalProblemDiff2d()
{
    // nothing here right now
    return; 
}

//---------------------------------------------------------------------------//
/*!
 * \brief This function returns an array of ResponseFunctionDiffusion objects.
 *
 * to be completed
 *
 */
ResponseFunction **LocalProblemDiff2d::getResponseFunctions( scalar keff )
{
    scalar t1 = MPI_Wtime();

    if ( keff == problem[0]->kNew ) // we're using the most recent
    {
        for (integer i = 0; i < input.numel; i++)
        {
            RFs[i] = new ResponseFunctionDiffusion( problem[i]->R, 
                                                    problem[i]->RL, 
                                                    problem[i]->RF, 
                                                    problem[i]->RA, 
                                                    problem[i]->kNew );
        }
        return RFs;
    }

    if (  keff == problem[0]->kOld ) // we're using the second most recent
    {
        for (integer i = 0; i < input.numel; i++)
        {
            RFs[i] = new ResponseFunctionDiffusion( problem[i]->Rold, 
                                                    problem[i]->RLold, 
                                                    problem[i]->RFold, 
                                                    problem[i]->RAold, 
                                                    problem[i]->kOld );
        }
        return RFs;
    }

    for (integer i = 0; i < input.numel; ++i)
    {

        if (input.ptype==2)
        {
            // initializes and/or resets the associated rf arrays
            problem[i]->resizeR( input.numg * (input.maxOrder+1) , keff );
            for (integer o = 0; o <= input.maxOrder; o++)    // incident orders
            {
                for (integer s = 0; s < 4; s++)              // incident faces
                {
                    for (integer g = 0; g < input.numg; g++) // incident groups
                    {
                        // produce the appropriate boundary source
                        problem[i]->updateRHS( keff, o, g, s );                 
                        // get the flux
                        solver[i]->solve();  
                        // compute the responses for g/o/s
                        solver[i]->compRespFct(); 
                        //solver[i]->destroy();
                    }
                }     
            }
        } // END BLOCK
        else
        {
            cout << " ERROR -- PTYPE =~ 2!! " << endl;
        }
    }
    for (integer i = 0; i < input.numel; i++)
    {
        RFs[i] = new ResponseFunctionDiffusion( problem[i]->R, 
                                                problem[i]->RL, 
                                                problem[i]->RF, 
                                                problem[i]->RA, 
                                                problem[i]->kNew );
    }
    keffOld = keffCurrent;
    keffCurrent = keff;

    scalar t2 = MPI_Wtime();
    time = time + (t2-t1);

    return RFs;
}
//---------------------------------------------------------------------------//
//                 end of LocalProblemDiff2d.cc
//---------------------------------------------------------------------------//

