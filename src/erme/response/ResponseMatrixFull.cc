//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseMatrixFull.cc
 * \author Jeremy Roberts
 * \date   11/24/2010
 * \brief  Member definitions of base class ResponseMatrixFull
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
#include "ResponseMatrix.hh"
#include "ResponseMatrixFull.hh"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
ResponseMatrixFull::ResponseMatrixFull( GlobalInput &input,
                                        ResponseFunctionServer *s )
 : ResponseMatrix(input,s), 
   SermentMatrixBCRS( numGroups*(1+spatialOrder)*(1+angularOrder)*numFaces*numElements,
                      numGroups*(1+spatialOrder)*(1+angularOrder)*numFaces*numElements, 
                      1, 
                      numGroups*(1+spatialOrder)*(1+angularOrder)*numFaces )
{
    std::cout << " CONSTRUCTING ResponseMatrixFull " << std::endl;
    elPlace(input);
    // nothing more here for now
}


//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
ResponseMatrixFull::~ResponseMatrixFull()
{
    // nothing more here right now, as M is destroyed explicitly in Release()
    return; 
}

//---------------------------------------------------------------------------//
/*!
 * \brief This function updates underlying matrix elements.
 *
 */
void ResponseMatrixFull::updateData( scalar k )
{
    // first update the response functions
    updateRF( k );
    // now, loop through element-wise and update the matrix blocks
    integer idx[1];
    integer rfIdx;
    integer numBlocks = 1;
    for ( integer i = 0; i < numElements; i++ )
    {
        idx[0] = i;
        rfIdx = elementPlacement[i];
        insertVals( responseFunctions[rfIdx]->getCurrentResponse(), 
                    numBlocks, idx, numBlocks, idx );
    }
    checkReady();
}

//---------------------------------------------------------------------------//
//                 end of ResponseMatrixFull.cc
//---------------------------------------------------------------------------//

