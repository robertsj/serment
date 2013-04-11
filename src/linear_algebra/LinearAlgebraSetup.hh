//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   LinearAlgebraSetup.hh
 *  @brief  Setup routines for linear algebra subsystems
 *  @author Jeremy Roberts
 *  @date   Sep 3, 2012
 */
//---------------------------------------------------------------------------//

#ifndef linear_algebra_LINEARALGEBRASETUP_HH_
#define linear_algebra_LINEARALGEBRASETUP_HH_

#include "comm/Comm.hh"
#include "petsc.h"
#include "slepc.h"

namespace linear_algebra
{

/// Initialize a parallel job.
void initialize(int &argc, char **&argv)
{
  // Require that comm is initialized
  Insist(serment_comm::Comm::is_comm_built(),
    "The local and global communicator must be built before linear algebra");

  // Set PETSc communicator to global if running in parallel
  if (serment_comm::Comm::is_global())
    PETSC_COMM_WORLD = serment_comm::global;
  else
    PETSC_COMM_WORLD = PETSC_COMM_SELF;

  // Initialize PETSc on the *global* communicator
  if (1)//serment_comm::Comm::is_global())
  {
    PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
    SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
  }

}

/// Finish a parallel job.
void finalize()
{
  // Finalize PETSc on the *global* communicator
  if (1)//(serment_comm::Comm::is_global())
  {
    SlepcFinalize();
    PetscFinalize();
  }
}


} // end namespace linear_algebra

#endif // linear_algebra_LINEARALGEBRASETUP_HH_

//---------------------------------------------------------------------------//
//              end of file LinearAlgebraSetup.hh
//---------------------------------------------------------------------------//
