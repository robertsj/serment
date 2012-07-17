//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseFunctionDiffusion.cc
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
#include <iostream>
#include "ResponseFunctionDiffusion.hh"
#include "linalg/typedefs.hh"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
ResponseFunctionDiffusion::ResponseFunctionDiffusion( scalar *C, scalar *L, 
                                                      scalar *F, scalar *A, scalar k )
  : ResponseFunction(C,L,F,A,k)
{
  //  std::cout << "hello responsefunctiondiffusion " << std::endl;
    // nothing more here for now
}

//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
ResponseFunctionDiffusion::~ResponseFunctionDiffusion()
{
    // nothing here right now
    return; 
}

//---------------------------------------------------------------------------//
//                 end of ResponseFunctionDiffusion.cc
//---------------------------------------------------------------------------//

