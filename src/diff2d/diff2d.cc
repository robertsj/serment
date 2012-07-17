//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   diff2d.cc
 * \author Jeremy Roberts
 * \date   08/26/2010
 * \brief  Driver for diff2d code
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 85                                            $:Rev of last commit
// $Author:: bert                                       $:Author of last commit
// $Date:: 2011-04-09 22:24:31 -0400 (Sat, 09 Apr 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#include <iostream>
#include <iomanip>
#include "Diff2dInput.hh"
#include "Diff2dProblem.hh"
#include "Diff2dSolver.hh"
#include "Diff2dOutput.hh"
#include "../linalg/typedefs.hh"
#include "petscmat.h"
using namespace std;
//---------------------------------------------------------------------------//

int main(int argc, char *args[])
{
    // initialize petsc
    PetscInitialize(&argc,&args,PETSC_NULL,PETSC_NULL);

    scalar t1, t2;

    t1 = MPI_Wtime();

    // create input
    Diff2dInput input;
    if ( argc < 2 ) cout << " error! need at least 1 command line argument. " << endl;
    input.readInput( args[1] );

    // create output
    Diff2dOutput output( input );

    // create the problem, setup and solve for each element
    Diff2dProblem *problem[ input.numel ];
    Diff2dSolver  *solver[ input.numel ];
    char tmp;
    for (integer i = 0; i < input.numel; ++i)
    {
        cout << " ================= PROBLEM " << i+1 << " ================= " << endl;
        problem[i] = new Diff2dProblem( input, i );

		// if this is a response function problem (for testing)
        if (input.ptype==2)
        {
            scalar keff = 1.0;
            problem[i]->resizeR( input.numg * (input.maxOrder+1) , 0.0);

            for (integer o = 0; o <= input.maxOrder; o++)     // for all orders
            {
                for (integer s = 0; s < 4; s++)               // for all faces
                {
                    for (integer g = 0; g < input.numg; g++)  // for all groups
                    {
                        problem[i]->updateRHS( keff, o, g, s );
                        solver[i]  = new Diff2dSolver( input, *problem[i], i );
                        // get the flux
                        solver[i]->solve();  
                        // compute the responses for g/o/s
                        solver[i]->compRespFct(); 
                        solver[i]->destroy();
						// Debug---watch the output
						//output.doOutput( *problem[i], *solver[i] );
						//output.closeSilo();
                        //cout << " enter any character " << endl;
                        //cin >> tmp;
                    }
                }     
            }

        } 
        else // stand-alone fixed-source or eigenvalue problem
        {
            solver[i]  = new Diff2dSolver( input, *problem[i], i );
            solver[i]->solve();            
            // do some output
            output.doOutput( *problem[i], *solver[i] );
        }
    }

    t2 = MPI_Wtime();
    cout << " solver time: " << t2 - t1 << " second " << endl;

    // kill things
    for (integer i = 0; i < input.numel; ++i)
    {
       // problem[i]->destroy();??
        delete problem[i];
        //solver[i]->destroy();
        delete solver[i];
    }

	PetscFinalize();
    return 0;
}

//---------------------------------------------------------------------------//
//                 end of diff2d.cc
//---------------------------------------------------------------------------//

