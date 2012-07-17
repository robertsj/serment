//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseFunctionDiffusion.hh
 * \author Jeremy Roberts
 * \date   11/24/2010
 * \brief  A container class for response function data.
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 140                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-09-14 12:53:40 -0400 (Wed, 14 Sep 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef RESPONSEFUNCTIONDIFFUSION_HH
#define RESPONSEFUNCTIONDIFFUSION_HH
#include "linalg/typedefs.hh"
#include "ResponseFunction.hh"

//===========================================================================//
/*!
 * \class ResponseFunctionDiffusion
 * \brief This base class is the foundation for holding and giving rf data.
 *
 * Each instance of ResponseFunctionDiffusion contains all the response 
 * function data for one element in a diffusion-based local problem.  To pass 
 * RF data of all elements, a pointer array of ResponseFunctionDiffusion 
 * objects is passed. Since we're using only pointers, this is pretty 
 * memory-efficient, I think.
 *
 */
//===========================================================================//

class ResponseFunctionDiffusion : public ResponseFunction
{

  public:
    // constructor
    ResponseFunctionDiffusion( scalar *C, scalar *L, scalar *F, scalar *A, scalar k );
    // destructor
    ~ResponseFunctionDiffusion();
    // Inherits:
    //    scalar *getCurrentResponse();
    //    scalar *getLeakageResponse();
    //    scalar *getFissionResponse();
    //    scalar *getAbsorptionResponse();
    //scalar keff;
  private:
    // Inherits:
    //    scalar *CurrentResponse;
    //    scalar *LeakageResponse;
    //    scalar *FissionResponse;
    //    scalar *AbsorptionResponse;
    
};

#endif // RESPONSEFUNCTIONDIFFUSION_HH

//---------------------------------------------------------------------------//
//                 end of ResponseFunctionDiffusion.hh
//---------------------------------------------------------------------------//

