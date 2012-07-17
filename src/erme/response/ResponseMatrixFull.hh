//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseMatrixFull.hh
 * \author Jeremy Roberts
 * \date   11/26/2010
 * \brief  A base class for response matrices.
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 167                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-11-07 23:01:04 -0500 (Mon, 07 Nov 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef RESPONSEMATRIXFULL_HH
#define RESPONSEMATRIXFULL_HH

#include "LinAlg.hh"
#include "ResponseFunction.hh"
#include "ResponseFunctionServer.hh"
#include "GlobalInput.hh"
#include "ResponseOperator.hh"
#include "ResponseMatrix.hh"

//===========================================================================//
/*!
 * \class ResponseMatrixFull
 * \brief This base class is part of the foundation of a whole response matrix.
 *
 * This class contains items specific to response matrices.  However, since we
 * wish to leverage the various matrix classes (full versus shell), we need 
 * multiple inheritance.
 *
 */
//===========================================================================//

class ResponseMatrixFull : public ResponseMatrix, public SermentMatrixBCRS
{

public:

  /// Typedefs
  //\{
  typedef typename util::SP<ResponseMatrixFull> SP_R;
  //\}

  // constructor
  ResponseMatrixFull(GlobalInput &input, ResponseFunctionServer *s);
  // destructor
  ~ResponseMatrixFull();
  // update the response matrix (needs specific matrix update)
  void updateData(scalar k);

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

#endif // ResponseMatrixFull_HH
//---------------------------------------------------------------------------//
//                 end of ResponseMatrixFull.hh
//---------------------------------------------------------------------------//

