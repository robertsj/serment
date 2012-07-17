//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseMatrix.cc
 * \author Jeremy Roberts
 * \date   11/24/2010
 * \brief  Member definitions of base class ResponseMatrix
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 140                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-09-14 12:53:40 -0400 (Wed, 14 Sep 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#include "LinAlg.hh"
#include "ResponseFunction.hh"
#include "ResponseFunctionServer.hh"
#include "GlobalInput.hh"
#include "ResponseOperator.hh"
#include "ResponseMatrix.hh"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
ResponseMatrix::ResponseMatrix( GlobalInput &input, 
                                ResponseFunctionServer *s )
  : ResponseOperator(input,s)
{
    // nothing more here for now
}

//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
ResponseMatrix::~ResponseMatrix()
{
    missServer = NULL;
    // nothing more here right now, as M is destroyed explicitly in Release()
    return; 
}

//---------------------------------------------------------------------------//
//                 end of ResponseMatrix.cc
//---------------------------------------------------------------------------//

