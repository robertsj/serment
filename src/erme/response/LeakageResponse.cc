//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LeakageResponse.cc
 * \author Jeremy Roberts
 * \date   11/24/2010
 * \brief  Member definitions of base class LeakageResponse
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
#include "LeakageResponse.hh"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//


LeakageResponse::LeakageResponse( GlobalInput &input,
                                  ResponseFunctionServer *s,
                                  integer *idx,              
                                  integer mmnum  )
 : ResponseMatrix(input,s), 
   SermentMatrixBCRS( numFaces*numElements,
                      numGroups*(1+spatialOrder)*(1+angularOrder)*numFaces*numElements, 
                      numGroups*(1+spatialOrder)*(1+angularOrder), 
                      numFaces ),
   mmIndex(idx), mmNum(mmnum)
{
    std::cout << " CONSTRUCTING LeakageResponse " << std::endl;
    elPlace(input);
    buildMM();
    // nothing more here for now
}


//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
LeakageResponse::~LeakageResponse()
{
    // nothing more here right now, as M is destroyed explicitly in Release()
    return; 
}

//---------------------------------------------------------------------------//
/*!
 * \brief This function updates underlying matrix elements.
 *
 */
void LeakageResponse::updateData( scalar k )
{
    // first update the response functions
    updateRF( k );
    // now, loop through element-wise and update the matrix blocks
    integer numBlocks = numGroups*(spatialOrder+1)*(angularOrder+1);
    integer idxy[numBlocks];
    integer idxx[1];
    integer rfIdx;
    integer one = 1;
    for ( integer i = 0; i < numElements; i++ )
    {
        for ( integer j = 0; j < numBlocks; j++ )
            idxy[j] = i*numBlocks+j;
        idxx[0] = i;
        rfIdx = elementPlacement[i];
        insertVals( responseFunctions[rfIdx]->getLeakageResponse(), 
                    one, idxx, numBlocks, idxy );
    }
    checkReady();
}

//---------------------------------------------------------------------------//
/*!
 * \brief This builds the zeroth-order connectivity operator.
 *
 */
void LeakageResponse::buildMM()
{
    
    if ( mmNum == 0 ) return; // no leakage, so don't bother making
 
    MM = new SermentVector( numFaces*numElements );
    result = new SermentVector( numFaces*numElements );

    scalar one = 1.0;
    for ( integer i = 0; i < mmNum; i++ )
    {
       // std::cout << "mmi= " << mmIndex[i] << std::endl;
        MM->insertVal( mmIndex[i], one );
    }
    //std::cout << " MM = " << std::endl;
    //MM->viewMe();
}

//---------------------------------------------------------------------------//
/*!
 * \brief This computes the leakage using the newest response function data.
 *
 */
scalar LeakageResponse::computeLeakage( SermentVector &J )
{
    if ( mmNum == 0 ) return 0.0; // no faces next to vacuum --> no leakage

    // result = L*J
    this->matVec( J, *result );
    // leakage = result'MM, ie result dot MM
    return MM->vecDot( *result );
}

void LeakageResponse::getLeakageVec( SermentVector &leakV )
{
    // net leakage is MM*(L*J), where MM is m*1, L is m*n, and J is n*1
    // we want the leakage vector with which to dot J.  This is simply
    // L'*MM which is n*m*m*1 = n*1, as needed
    MatMultTranspose( M, MM->V, leakV.V );
    return;
}


//---------------------------------------------------------------------------//
//                 end of LeakageResponse.cc
//---------------------------------------------------------------------------//

