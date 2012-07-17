//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   AbsorptionResponse.hh
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

#ifndef ABSORPTIONRESPONSE_HH
#define ABSORPTIONRESPONSE_HH
#include "LinAlg.hh"
#include "ResponseFunction.hh"
#include "ResponseFunctionServer.hh"
#include "GlobalInput.hh"
#include "ResponseOperator.hh"

//===========================================================================//
/*!
 * \class AbsorptionResponse
 * \brief This class encapsulates a SermentVector and required RF updates.
 *
 */
//===========================================================================//

class AbsorptionResponse : public ResponseOperator, public SermentVector
{

  public:
    // constructor
    AbsorptionResponse( GlobalInput &input,
                        ResponseFunctionServer *s );
    // destructor
    ~AbsorptionResponse();
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

#endif // ABSORPTIONRESPONSE_HH

//---------------------------------------------------------------------------//
//                 end of AbsorptionResponse.hh
//---------------------------------------------------------------------------//

