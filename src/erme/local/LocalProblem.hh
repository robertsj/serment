//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LocalProblem.hh
 * \author Jeremy Roberts
 * \date   11/23/2010
 * \brief  Base class for local problem routines.
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 140                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-09-14 12:53:40 -0400 (Wed, 14 Sep 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef LOCALPROBLEM_HH
#define LOCALPROBLEM_HH

#include <iostream>
#include "ResponseFunction.hh"
#include "linalg/typedefs.hh"

using namespace std;
//===========================================================================//
/*!
 * \class LocalProblem
 * \brief This class encapsulates the entire local problem routine.
 *
 * to be completed
 *
 */
//===========================================================================//

class LocalProblem
{
  public:

    LocalProblem( string in ) : LocalInputFile(in) { time=0; };
    ~LocalProblem(){};
    virtual ResponseFunction **getResponseFunctions( scalar keff ) = 0;
    scalar localTime(){ return time; };

  private:
    string LocalInputFile;

  protected:
    scalar time;


};

#endif // LOCALPROBLEM_HH

//---------------------------------------------------------------------------//
//                 end of LocalProblem.hh
//---------------------------------------------------------------------------//

