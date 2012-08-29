//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Serial.cc
 * \brief  Serial comm implementation member definitions.
 * \author Jeremy Roberts
 * \date   Aug 21, 2012
 */
//---------------------------------------------------------------------------//

#include "Comm.hh"

#ifndef SERMENT_ENABLE_MPI

namespace serment_comm
{

// Initialize private static variables
bool   Comm::d_global = false;
bool   Comm::d_local  = false;
double Comm::d_time   = 0.0;

//---------------------------------------------------------------------------//
// SETUP FUNCTIONS
//---------------------------------------------------------------------------//

void Comm::initialize(int &argc, char **&argv)
{
  /* ... */
}

void Comm::finalize()
{
  /* ... */
}

//---------------------------------------------------------------------------//
// COMMUNICATORS
//---------------------------------------------------------------------------//

/*!
 *  \brief Create communicators
 *  \param N  Number of local groups to create
 */
void Comm::setup_communicators(const unsigned int N)
{
  /* ... */
}

static bool Comm::is_global = true;

//---------------------------------------------------------------------------//
// QUERY FUNCTIONS
//---------------------------------------------------------------------------//

int rank()
{
  return 0;
}

int size()
{
  return 1;
}

} // end namespace serment_comm

#endif

//---------------------------------------------------------------------------//
//              end of file Serial.cc
//---------------------------------------------------------------------------//
