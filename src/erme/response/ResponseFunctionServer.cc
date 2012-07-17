//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseFunctionServer.cc
 * \author Jeremy Roberts
 * \date   11/24/2010
 * \brief  Member definitions of base class ResponseFunctionServer
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 140                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-09-14 12:53:40 -0400 (Wed, 14 Sep 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#include "iostream"
#include "LinAlg.hh"
#include "GlobalInput.hh"
#include "ResponseFunction.hh"
#include "DataBase.hh"
#include "LocalProblem.hh"
#include "LocalProblemDiff2d.hh"
#include "ResponseFunctionServer.hh"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
ResponseFunctionServer::ResponseFunctionServer( GlobalInput &input )
{
    currentKeff = -100.0;
    oldKeff = -200.0;
    time = 0.0;

    // INITIATE RF SOURCE
    if ( input.rfsource == 1 )
    {
        localGuy = new LocalProblemDiff2d( input.rfsourcefile );
    }
    else
    {
        std::cout << " ERROR: rfsource=" << input.rfsource 
                  << " NOT impemented " << std::endl;
    }
    // INITIATE RF SOURCE
    if ( input.rfsink == 1 )
    {
        std::cout << " PERFORMING GLOBAL PROBLEM " << std::endl;
        // localGuy = new LocalProblemDiff2d;
    }
    else
    {
        std::cout << " ERROR: rfsink=" << input.rfsource 
                  << " NOT impemented " << std::endl;
    }


}

//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
ResponseFunctionServer::~ResponseFunctionServer()
{
    // nothing more here right now
}

//---------------------------------------------------------------------------//
/*!
 * \brief This function updates the response functions via the server.
 *
 */
ResponseFunction **ResponseFunctionServer::updateResponseFunctions( scalar k )
{
        scalar t1 = MPI_Wtime(); 

        // eventually, there will be more here when a database class is
        // implemented; for now, rf's come only from local solver
        currentResponseFunction = localGuy->getResponseFunctions( k );

        scalar t2 = MPI_Wtime();
        time = time + (t2-t1);

        return currentResponseFunction;

}

//---------------------------------------------------------------------------//
//                 end of ResponseFunctionServer.cc
//---------------------------------------------------------------------------//

