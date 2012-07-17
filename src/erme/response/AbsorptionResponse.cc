//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   AbsorptionResponse.cc
 * \author Jeremy Roberts
 * \date   11/27/2010
 * \brief  Member definitions of base class AbsorptionResponse
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 140                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-09-14 12:53:40 -0400 (Wed, 14 Sep 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#include <iostream>
#include "LinAlg.hh"
#include "ResponseFunction.hh"
#include "ResponseFunctionServer.hh"
#include "GlobalInput.hh"
#include "ResponseOperator.hh"
#include "AbsorptionResponse.hh"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
AbsorptionResponse::AbsorptionResponse( GlobalInput &input,   
                                        ResponseFunctionServer *s )
 : ResponseOperator(input,s), 
   SermentVector(numGroups*(1+spatialOrder)*(1+angularOrder)*numFaces*numElements)
{
    std::cout << " CONSTRUCTING AbsorptionResponse " << std::endl;
    elPlace(input);
    // nothing more here for now
}

//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
AbsorptionResponse::~AbsorptionResponse()
{
    // nothing more here right now, as M is destroyed explicitly in Release()
    return; 
}

//---------------------------------------------------------------------------//
/*!
 * \brief This function updates underlying matrix elements.
 *
 */
void AbsorptionResponse::updateData( scalar k )
{
    // first update the response functions
    responseFunctions = missServer->updateResponseFunctions( k ); 
    // now, loop through element-wise and update the matrix blocks
    integer len = numGroups*(1+spatialOrder)*(1+angularOrder)*4;
    integer idx[len];
    integer rfIdx;
    scalar *abs;
//void SermentVector::insertVals( integer ni, const integer ix[],  
//                                const scalar y[] )
    for ( integer i = 0; i < numElements; i++ )
    {
        for ( integer j = 0; j < len; j++ )
            idx[j] = i*len+j; 
        rfIdx = elementPlacement[i];
        abs = responseFunctions[rfIdx]->getAbsorptionResponse();
        insertVals( len, idx, abs );
    }
    checkReady();
}

//---------------------------------------------------------------------------//
//                 end of AbsorptionResponse.cc
//---------------------------------------------------------------------------//

