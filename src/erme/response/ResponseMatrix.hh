//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseMatrix.hh
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

#ifndef RESPONSEMATRIX_HH
#define RESPONSEMATRIX_HH

#include "LinAlg.hh"
#include "ResponseFunction.hh"
#include "ResponseFunctionServer.hh"
#include "GlobalInput.hh"
#include "ResponseOperator.hh"

//===========================================================================//
/*!
 * \class ResponseMatrix
 * \brief This base class is part of the foundation of a whole response matrix.
 *
 * This class contains items specific to response matrices.  However, since we
 * wish to leverage the various matrix classes (full versus shell), we need 
 * multiple inheritance.
 *
 */
//===========================================================================//

class ResponseMatrix : public ResponseOperator
{

  public:
    // constructor
    ResponseMatrix( GlobalInput &input,
                    ResponseFunctionServer *s );
    // destructor
    ~ResponseMatrix();
    // update the response matrix (needs specific matrix update)
    //virtual void updateData( scalar k ) = 0;

  protected:
    // inherits
    //   integer numGroups;
    //   integer spatialOrder;
    //   integer angularOrder;
    //   integer numFaces;
    //   integer numElements;
    //   integer *elementPlacement;
    //   ResponseFunctionServer *missServer;
    //   ResponseFunction **responseFunctions;

    
};

#endif // RESPONSEMATRIX_HH

//---------------------------------------------------------------------------//
//                 end of ResponseMatrix.hh
//---------------------------------------------------------------------------//

