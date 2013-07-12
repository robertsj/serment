//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  xerme.cc
 *  @brief Serment executable driver
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "serment_config.h"
#include "erme_solver/ManagerERME.hh"
#include "erme_utils/Archive.hh"
#include "utilities/Timer.hh"
#include <iostream>
#include <ctime>

void print_welcome();

//----------------------------------------------------------------------------//
int main(int argc, char **argv)
{
  Insist(argc > 1, "Not enough command line arguments!");

  // Setup the communicators
  erme_solver::ManagerERME M(argc, argv);

  print_welcome();

  // Get input from archive.  This assumes the first argument is
  // the archive file name.
  erme_utils::Archive A;
  A.unarchive(std::string(argv[1]));

  // Setup and solve problem
  M.build_comm(A.get_db());
  M.build_erme(A.get_nodes());
  M.solve();

  return 0;
}

//----------------------------------------------------------------------------//
void print_welcome()
{
  using std::cout;
  using std::endl;
  if (serment_comm::Comm::world_rank() == 0)
  {
    std::time_t t;
    std::time(&t);
    cout << "                                    _       " << endl;
    cout << "                                   | |      " << endl;
    cout << " ___  ___ _ __ _ __ ___   ___ _ __ | |_     " << endl;
    cout << "/ __|/ _ \\ '__| '_ ` _ \\ / _ \\ '_ \\| __|" << endl;
    cout << "\\__ \\  __/ |  | | | | | |  __/ | | | |_   " << endl;
    cout << "|___/\\___|_|  |_| |_| |_|\\___|_| |_|\\__| " << endl;
    cout <<
      " solving eigenvalue response matrix equations w/ nonlinear techniques"
         << endl;
    cout << " Copyright (C) 2012-2013 Jeremy Roberts"  << endl << endl;
    cout << " Serment built on: " << SERMENT_COMPILED_M << "/"
                                  << SERMENT_COMPILED_D << "/"
                                  << SERMENT_COMPILED_Y << endl;
    cout << "         Git SHA1: " << SERMENT_GIT_SHA1  << endl;
    cout << "  Detran built on: " << DETRAN_COMPILED_M << "/"
                                  << DETRAN_COMPILED_D << "/"
                                  << DETRAN_COMPILED_Y << endl;
    cout << "         Git SHA1: " << DETRAN_GIT_SHA1  << endl;
    cout << "           Run on: " << std::ctime(&t) << endl << endl;
  }
}

//----------------------------------------------------------------------------//
//              end of file xerme.cc
//----------------------------------------------------------------------------//
