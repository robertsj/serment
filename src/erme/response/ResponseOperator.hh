//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseOperator.hh
 * \author Jeremy Roberts
 * \date   11/26/2010
 * \brief  A base class for response matrices.
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 140                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-09-14 12:53:40 -0400 (Wed, 14 Sep 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef RESPONSEOPERATOR_HH
#define RESPONSEOPERATOR_HH

#include <iostream>
#include "LinAlg.hh"
#include "ResponseFunction.hh"
#include "ResponseFunctionServer.hh"
#include "GlobalInput.hh"

//===========================================================================//
/*!
 * \class ResponseOperator
 * \brief This base class has general traits of a response operator.
 *
 * Such operators include the response matrix, leakage response matric, etc.
 *
 */
//===========================================================================//

class ResponseOperator
{

  public:
    // constructor
    ResponseOperator( GlobalInput &input,
                      ResponseFunctionServer *s )
        : numGroups(input.numgroups), spatialOrder(input.spaceord),
          angularOrder(input.angleord), numFaces(input.faces),
          numElements(input.numel), missServer(s) {};
    // destructor
    ~ResponseOperator(){};

    // update the underlying response data for the given eigenvalue
    virtual void updateData( scalar k ) = 0;

  protected:
    /*!
     * \brief number of groups
     */
    integer numGroups;
    /*!
     * \brief spatial expansion order
     */
    integer spatialOrder;
    /*!
     * \brief angular expansion order
     */
    integer angularOrder;
    /*!
     * \brief number of element faces
     */
    integer numFaces;
    /*!
     * \brief number of elements (coarse meshes) comprising global problem
     */
    integer numElements;
    /*!
     * \brief placement of the various element types ( ordered by... ?)
     */
    integer *elementPlacement;
    /*!
     * \brief pointer to the server to whom requests are made for rf-updates
     */
    ResponseFunctionServer *missServer;
    /*!
     * \brief pointers to constituent rf's
     */
    ResponseFunction **responseFunctions;
    /*!
     * \brief This function updates the response functions via the server.
     *
     */
    void updateRF( scalar k )
    {
        responseFunctions = missServer->updateResponseFunctions( k );    
    }
    /*!
     * \brief This builds the elementPlacement from input.
     *
     * Note, this counts elements by physical row, beginning with
     * the bottom, for a single z-place.  That means we go along
     * x, then increment y, and repeat for all z.
     */
    void elPlace(GlobalInput &input)
    {
        elementPlacement = new integer[input.numel];
        integer e = 0;
        for ( integer k=0; k < input.elemz; k++)
        {
            for ( integer j=0; j < input.elemy; j++)
            {
                for ( integer i=0; i < input.elemx; i++)
                {
                    if ( input.elements[i][j][k] >= 0 )
                    {
                        elementPlacement[e] = input.elements[i][j][k];
                        e++;
                    }
                }
            }
        }   
        if ( e != input.numel )
            std::cout << " e=" << e << " ~= numel=" << input.numel << std::endl; 

    }
};

#endif // RESPONSEOPERATOR_HH

//---------------------------------------------------------------------------//
//                 end of ResponseOperator.hh
//---------------------------------------------------------------------------//

