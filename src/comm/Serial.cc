//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Serial.cc
 *  @brief Serial comm implementation member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "Comm.hh"

#ifndef SERMENT_ENABLE_MPI

namespace serment_comm
{

// Set the default communicators
Communicator_t world          = 0;
Communicator_t global         = world;
Communicator_t local          = world;
Communicator_t communicator   = world;

// Initialize private static variables
bool   Comm::d_is_comm_built = false;
bool   Comm::d_is_global     = true;
double Comm::d_time = 0.0;
int    Comm::d_world_rank = 0;
int    Comm::d_local_group = 0;

//----------------------------------------------------------------------------//
// SETUP FUNCTIONS
//----------------------------------------------------------------------------//

void Comm::initialize(int &argc, char **&argv)
{
  /* ... */
}

void Comm::finalize()
{
  /* ... */
}

//----------------------------------------------------------------------------//
// COMMUNICATORS
//----------------------------------------------------------------------------//

/**
 *  @brief Create communicators
 *  @param N  Number of local groups to create
 */
void Comm::setup_communicators(const unsigned int N)
{
  d_is_comm_built = true;
}

//----------------------------------------------------------------------------//
// QUERY FUNCTIONS
//----------------------------------------------------------------------------//

int Comm::rank()
{
  return 0;
}

int Comm::size()
{
  return 1;
}

} // end namespace serment_comm

#endif

//----------------------------------------------------------------------------//
//              end of file Serial.cc
//----------------------------------------------------------------------------//
