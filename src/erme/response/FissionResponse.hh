//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   FissionResponse.hh
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

#ifndef FISSIONRESPONSE_HH
#define FISSIONRESPONSE_HH
#include "LinAlg.hh"
#include "ResponseFunction.hh"
#include "ResponseFunctionServer.hh"
#include "GlobalInput.hh"
#include "ResponseOperator.hh"

//===========================================================================//
/*!
 * \class FissionResponse
 * \brief This class encapsulates a SermentVector and required RF updates.
 *
 */
//===========================================================================//

class FissionResponse : public ResponseOperator, public SermentVector
{

  public:
    // constructor
    FissionResponse( GlobalInput &input,
                     ResponseFunctionServer *s );
    // destructor
    ~FissionResponse();
    // update the response matrix (needs specific matrix update)
    void updateData( scalar k );
  private:
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

#endif // FISSIONRESPONSE_HH

//---------------------------------------------------------------------------//
//                 end of FissionResponse.hh
//---------------------------------------------------------------------------//

