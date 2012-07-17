//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LeakageResponse.hh
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

#ifndef LEAKAGERESPONSE_HH
#define LEAKAGERESPONSE_HH
#include "LinAlg.hh"
#include "ResponseFunction.hh"
#include "ResponseFunctionServer.hh"
#include "GlobalInput.hh"
#include "ResponseOperator.hh"
#include "ResponseMatrix.hh"

//===========================================================================//
/*!
 * \class LeakageResponse
 * \brief This is a matrix operator for computing the leakage response.
 *
 * Because it is a matrix quantity similar to the response matrix, we build
 * from the same common base and use a BCRS matrix format.  Note, the block
 * size is the number of edges/faces.  For example, a 2d diffusion element
 * has a leakage response function comprised of four vectors, one for each 
 * side.  The inner product of these vectors with the incident current vector
 * (incident on that element) produce the net leakage from the element.
 *
 */
//===========================================================================//

class LeakageResponse : public ResponseMatrix, public SermentMatrixBCRS
{

  public:
    // constructor
    LeakageResponse( GlobalInput &input,
                     ResponseFunctionServer *s,
                     integer *idx,             
                     integer mmnum  );
    // destructor
    ~LeakageResponse();
    // update the response matrix (needs specific matrix update)
    void updateData( scalar k );
    scalar computeLeakage( SermentVector &J );
    void getLeakageVec( SermentVector &leakV );
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
    integer *mmIndex;
    integer mmNum;
    SermentVector *MM; // "little M", ones exist where leakage occurs; dot with L*J+ to get total leakage
    SermentVector *result;
    void buildMM();

};

#endif // LeakageResponse_HH

//---------------------------------------------------------------------------//
//                 end of LeakageResponse.hh
//---------------------------------------------------------------------------//

