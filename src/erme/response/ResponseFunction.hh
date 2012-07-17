//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseFunction.hh
 * \author Jeremy Roberts
 * \date   11/24/2010
 * \brief  A container class for response function data.
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 157                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-10-14 17:05:41 -0400 (Fri, 14 Oct 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef RESPONSEFUNCTION_HH
#define RESPONSEFUNCTION_HH
#include "LinAlg.hh"

//===========================================================================//
/*!
 * \class ResponseFunction
 * \brief This base class is the foundation for holding and giving rf data.
 *
 * Each instance of ResponseFunction (or more specifically, it's subclasses)
 * contains all the response function data for one element.  To pass RF data
 * of all elements, a pointer array of ResponseFunction objects is passed.
 * Since we're using only pointers, this is pretty memory-efficient, I think.
 *
 */
//===========================================================================//

class ResponseFunction
{

public:
  // constructor
  ResponseFunction(scalar *C, scalar *L, scalar *F, scalar *A, scalar k);
  // destructor
  ~ResponseFunction();
  // We keep the RF's private, and return the associated pointer with these:
  scalar *getCurrentResponse();
  scalar *getLeakageResponse();
  scalar *getFissionResponse();
  scalar *getAbsorptionResponse();

  // private:
  // Baseline data -- needed for all potential responses.
  // We use pointers so that no copying of data is needed, i.e. we use the
  //   physical data produced via the LocalProblem.
  scalar *CurrentResponse;
  scalar *LeakageResponse;
  scalar *FissionResponse;
  scalar *AbsorptionResponse;
  scalar keff;
};

#endif // RESPONSEFUNCTION_HH

//---------------------------------------------------------------------------//
//                 end of ResponseFunction.hh
//---------------------------------------------------------------------------//

