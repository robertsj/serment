//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LinearAlgebraSetup.hh
 * \brief  Setup routines for linear algebra subsystems
 * \author Jeremy Roberts
 * \date   Sep 3, 2012
 */
//---------------------------------------------------------------------------//

#ifndef LINEARALGEBRASETUP_HH_
#define LINEARALGEBRASETUP_HH_

#include "comm/Comm.hh"
#include "petsc.h"

namespace linear_algebra
{

/// Initialize a parallel job.
void initialize(int &argc, char **&argv)
{
  // Require that comm is initialized
  Insist(serment_comm::Comm::is_comm_built(),
    "The local and global communicator must be built before linear algebra");

  // Set PETSc communicator to global
  PETSC_COMM_WORLD = serment_comm::global;

  // Initialize PETSc
  PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
}

/// Finish a parallel job.
void finalize()
{
  PetscFinalize();
}


} // end namespace linear_algebra

#endif // LINEARALGEBRASETUP_HH_ 

//---------------------------------------------------------------------------//
//              end of file LinearAlgebraSetup.hh
//---------------------------------------------------------------------------//
