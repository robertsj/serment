//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Serment.cc
 * \author Jeremy Roberts
 * \date   11/05/2010
 * \brief  Executable driver for SERMENT.
 * \note   Copyright (C) 2011 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 178                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-12-12 13:36:39 -0500 (Mon, 12 Dec 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#include <iostream>
#include <cmath>

#include "serment_config.h"

#include "LinAlg.hh"
#include "GlobalInput.hh"
#include "ResponseFunctionServer.hh"
#include "GlobalProblem.hh"
#include "PowerIter.hh"
#include "Newton.hh"

// Inner iterations
#include "InnerIterPower.hh"
#ifdef SERMENT_ENABLE_SLEPC
#include "InnerIterSLEPc.hh"
#endif

#include "utilities/DBC.hh"

using namespace std;

//---------------------------------------------------------------------------//

int main(int argc, char *args[])
{

  Insist(argc > 1, "Need at least 1 command line argument.");

  PetscInitialize(&argc, &args, PETSC_NULL, PETSC_NULL);
#ifdef SERMENT_ENABLE_SLEPC
  SlepcInitialize(&argc, &args, (char*) 0, "");
#endif

  scalar tTOT = MPI_Wtime();

  // Global input with the grammar as an argument.
  GlobalInput::SP_globalinput input;
  input = new GlobalInput("serment.rng");

  if (!input->readInput(args[1]))
  {
    cout << " input error...goodbye " << endl;
    return 0;
  }
  input->echoInput();

  // Create the server; she decides to/from whom rf's go/come
  ResponseFunctionServer *s = new ResponseFunctionServer(*input);

  GlobalProblem::SP_globalproblem problem;
  problem = new GlobalProblem(input, s);

  scalar tPI = 0.0; // PI time
  scalar tN = 0.0; // Newton time

  // PI solve
  tPI = MPI_Wtime();
  if (input->pctype == -1)
  {
    PowerIter<InnerIterPower> solver(problem, input);
    solver.solve();
  }
  if (input->pctype == -2)
  {
#ifdef SERMENT_ENABLE_SLEPC
    PowerIter<InnerIterSLEPc> solver(problem, input);
    for (int i = 0; i < 100; i++)
        solver.solve();
#endif
  }
  // PI time.
  tPI = MPI_Wtime() - tPI;

  // Newton solve
  tN = MPI_Wtime();
  if (input->pctype >= 0)
  {
    Newton solver2(problem, input);
    for (int i = 0; i < 100; i++)
        solver2.solve();
  }

  // Compute Newton time
  tN = MPI_Wtime() - tN;

  // Compute total time
  tTOT = MPI_Wtime() - tTOT;

  std::cout << "  PI     Time: " << tPI << std::endl;
  std::cout << "  Newton Time: " << tN << std::endl;
  std::cout << "  Server Time: " << s->serverTime() << std::endl;
  std::cout << "  Total  Time: " << tTOT << std::endl;

  PetscFinalize();

  std::cout << "  ... te ootte nyt menneet loppuun " << std::endl;

  return 0;
}

//---------------------------------------------------------------------------//
//                 end of Serment.cc
//---------------------------------------------------------------------------//

