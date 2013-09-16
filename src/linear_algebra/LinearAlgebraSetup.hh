//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  LinearAlgebraSetup.hh
 *  @brief Setup routines for linear algebra subsystems
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef linear_algebra_LINEARALGEBRASETUP_HH_
#define linear_algebra_LINEARALGEBRASETUP_HH_

#include "comm/Comm.hh"
#include "petsc.h"
#include "slepc.h"

namespace linear_algebra
{

/// Initialize a parallel job.
void initialize(int &argc, char **&argv, bool init_comm = false)
{
  if (init_comm)
  {
    // Initialize the comm so that world = global, i.e. each process
    // is a local group.  This makes the most sense for testing.
    serment_comm::Comm::initialize(argc, argv);
    int N = serment_comm::Comm::size();
    serment_comm::Comm::setup_communicators(N);
  }

  Insist(serment_comm::Comm::is_comm_built(),
         "Communicators must be built before linear algebra");

  // Set PETSc communicator to global if running in parallel
#ifdef SERMENT_ENABLE_MPI
  if (serment_comm::Comm::is_global())
    PETSC_COMM_WORLD = serment_comm::global;
  else
    PETSC_COMM_WORLD = PETSC_COMM_SELF;
#else
  //PETSC_COMM_WORLD = PETSC_COMM_SELF;
#endif

  // Initialize PETSc on the *global* communicator
  if (1)//serment_comm::Comm::is_global())
  {
    PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
    SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
  }

}

/// Finish a parallel job.
void finalize(bool final_comm = false)
{
  // Finalize PETSc on the *global* communicator
  if (1)//(serment_comm::Comm::is_global())
  {
    SlepcFinalize();
    PetscFinalize();
  }

  if (final_comm)
  {
    serment_comm::Comm::finalize();
  }
}

} // end namespace linear_algebra

#endif // linear_algebra_LINEARALGEBRASETUP_HH_

//----------------------------------------------------------------------------//
//              end of file LinearAlgebraSetup.hh
//----------------------------------------------------------------------------//
