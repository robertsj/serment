//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseFunctionServer.hh
 * \author Jeremy Roberts
 * \date   11/26/2010
 * \brief  Routes response functions.
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 140                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-09-14 12:53:40 -0400 (Wed, 14 Sep 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef RESPONSEFUNCTIONSERVER_HH
#define RESPONSEFUNCTIONSERVER_HH
#include "LinAlg.hh"
#include "GlobalInput.hh"
#include "ResponseFunction.hh"
#include "DataBase.hh"
#include "LocalProblem.hh"


//===========================================================================//
/*!
 * \class ResponseFunctionServer
 * \brief To be completed
 *
 */
//===========================================================================//

class ResponseFunctionServer
{

  public:
    // constructor
    ResponseFunctionServer( GlobalInput &input );
    // destructor
    ~ResponseFunctionServer();
    // getResponseFunctions
    ResponseFunction **updateResponseFunctions( scalar k );
    // get server time
    scalar serverTime(){ return time; };
    // reset server time
    scalar resetTime(){ time = 0; };
  private:
    LocalProblem *localGuy;
    DataBase *dataGuy;
    scalar currentKeff;
    scalar oldKeff;
    ResponseFunction **currentResponseFunction;
    ResponseFunction **lastResponseFunction;
    ResponseFunction **tmp;
    scalar time;
};

#endif // RESPONSEFUNCTIONSERVER_HH

//---------------------------------------------------------------------------//
//                 end of ResponseFunctionServer.hh
//---------------------------------------------------------------------------//

