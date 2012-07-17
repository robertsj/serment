//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   GlobalSolver.cc
 * \author Jeremy Roberts
 * \date   11/24/2010
 * \brief  Member definitions of abstract class GlobalSolver
 * \note   Copyright (C) 2011 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 167                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-11-07 23:01:04 -0500 (Mon, 07 Nov 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#include <iostream>
#include <cmath>
#include "LinAlg.hh"
#include "GlobalInput.hh"
#include "ResponseFunctionServer.hh"
#include "ResponseMatrix.hh"
#include "ResponseMatrixFull.hh"
#include "AbsorptionResponse.hh"
#include "FissionResponse.hh"
#include "LeakageResponse.hh"
#include "ConnectMatrix.hh"
#include "Connect2dCart.hh"
#include "GlobalProblem.hh"
#include "GlobalSolver.hh"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
GlobalSolver::GlobalSolver( SP_globalproblem problem, SP_globalinput input )
 : d_problem(problem), d_input(input)
{
	return;
}

//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
GlobalSolver::~GlobalSolver()
{
    // nothing more here right now, as M is destroyed explicitly in Release()
    return; 
}


//---------------------------------------------------------------------------//
/*!
 * \brief This function computes and prints the coarse mesh fission rates.
 *
 * \param J    latest current coefficient vector
 */
void GlobalSolver::fissionRates( SP_vector J )
{

	// create a new vector for the fission rate
    SermentVector fissionrate( J->Length() );
	
	// compute the fission rate using point-wise multiplication
	J->vecPointMult( d_problem->F, fissionrate );

	// degrees of free per element
	integer dofperel = d_input->numgroups *
                       (1+d_input->spaceord)*(1+d_input->angleord)*d_input->faces;

	scalar *fissionrate_a;                   // fissionrate underlying array
	scalar elementfissionrate[d_input->numel]; // fission rate in each element
	int idx; // dummy index

	for ( int e = 0; e < d_input->numel; e++)
		elementfissionrate[e] = 0.0;

	VecGetArray( fissionrate.V, &fissionrate_a );  // get the array
	for ( int e = 0; e < d_input->numel; e++ )
	{
		for ( int i = 0; i < dofperel; i++ )
		{
			idx = e*dofperel + i;
			elementfissionrate[e] = elementfissionrate[e] + 
								    fissionrate_a[idx];
		}
	}
	VecRestoreArray( fissionrate.V, &fissionrate_a );  // must put array back

	idx = 0;
	for ( int k = 0; k < d_input->elemz; k++ )           // for all z-planes
	{
		for ( int i = 0; i < d_input->elemx; i++ )       // for all x
		{
			std::cout << std::endl;
			for ( int j = 0; j < d_input->elemy; j++ )    // for all y
			{
				if ( d_input->elements[i][j][k] == -1 )
					printf (" %8.6f", 0.0 );
				else
	 				printf (" %8.6f", elementfissionrate[idx] );
				idx++;
			}
		}
		std::cout << std::endl;	
	}
	// cleanup our mess
	fissionrate.releaseMe();
}


//---------------------------------------------------------------------------//
//                 end of GlobalSolver.cc
//---------------------------------------------------------------------------//

