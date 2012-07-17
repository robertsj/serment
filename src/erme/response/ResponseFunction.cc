//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseFunction.cc
 * \author Jeremy Roberts
 * \date   11/24/2010
 * \brief  Member definitions of base class ResposeFunction
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 140                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-09-14 12:53:40 -0400 (Wed, 14 Sep 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#include "ResponseFunction.hh"
#include "linalg/typedefs.hh"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
ResponseFunction::ResponseFunction( scalar *C, scalar *L, scalar *F, scalar *A, scalar k )
  : CurrentResponse(C), LeakageResponse(L), 
    FissionResponse(F), AbsorptionResponse(A), keff(k)
{
    // nothing more here for now
}

//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
ResponseFunction::~ResponseFunction()
{
    CurrentResponse = NULL;
    LeakageResponse = NULL;
    FissionResponse = NULL;
    AbsorptionResponse = NULL;
    // nothing more here right now, as M is destroyed explicitly in Release()
    return; 
}

//---------------------------------------------------------------------------//
/*!
 * \brief This function returns a pointer to current response data.
 *
 */
scalar *ResponseFunction::getCurrentResponse()
{
    return CurrentResponse;    
}

//---------------------------------------------------------------------------//
/*!
 * \brief This function returns a pointer to leakage response data.
 *
 */
scalar *ResponseFunction::getLeakageResponse()
{
    return LeakageResponse;    
}

//---------------------------------------------------------------------------//
/*!
 * \brief This function returns a pointer to fission response data.
 *
 */
scalar *ResponseFunction::getFissionResponse()
{
    return FissionResponse;    
}

//---------------------------------------------------------------------------//
/*!
 * \brief This function returns a pointer to absorption response data.
 *
 */
scalar *ResponseFunction::getAbsorptionResponse()
{
    return AbsorptionResponse;    
}


//---------------------------------------------------------------------------//
//                 end of ResponseFunction.cc
//---------------------------------------------------------------------------//

