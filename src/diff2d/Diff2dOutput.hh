//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Diff2dOutput.hh
 * \author Jeremy Roberts
 * \date   10/26/2010
 * \brief  A class for handling input for diff2d problems.
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 177                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-12-09 16:31:49 -0500 (Fri, 09 Dec 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef DIFF2DOUTPUT_HH
#define DIFF2DOUTPUT_HH

#include "serment_config.h"

#include "../linalg/typedefs.hh"

#include "Diff2dInput.hh"
#include "Diff2dProblem.hh"
#include "Diff2dSolver.hh"

#ifdef SERMENT_ENABLE_SILO
#include "silo.h"
#endif

using namespace std;

//===========================================================================//
/*!
 * \class Diff2dOutput
 * \brief A class for outputing fluxes and other things.
 *
 *  Diff2dOutput is blah blah blah.
 *
 */
//===========================================================================//

class Diff2dOutput
{
  public:
    
    Diff2dOutput( Diff2dInput &inp );
    ~Diff2dOutput();
    void doOutput( Diff2dProblem &prob, Diff2dSolver &sol );
    void closeSilo();

  private:
    bool silocreated;     // if the silo file is already created, don't repeat
    Diff2dInput *inp;     // pointer to the problem input
#ifdef SERMENT_ENABLE_SILO
    DBfile *file;         // The Silo file pointer   
#endif
    void printFlux();
    void plotFlux( Diff2dProblem &prob, Diff2dSolver &sol ); 

};

#endif // DIFF2DOUTPUT_HH

//---------------------------------------------------------------------------//
//                 end of Diff2dOutput.hh
//---------------------------------------------------------------------------//

