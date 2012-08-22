//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MPI.cc
 * \brief  MPI 
 * \author Jeremy Roberts
 * \date   Aug 21, 2012
 */
//---------------------------------------------------------------------------//

#include "Comm.hh"

namespace serment_comm
{

MPI_Comm default_communicator = MPI_COMM_WORLD;
MPI_Comm communicator = default_communicator;

//---------------------------------------------------------------------------//
// SETUP FUNCTIONS
//---------------------------------------------------------------------------//

double Comm::d_time = 0.0;

void Comm::initialize(int &argc, char **&argv)
{
  int result = MPI_Init(&argc, &argv);
  d_time = 0.0;
  Ensure(result == MPI_SUCCESS);
}

//---------------------------------------------------------------------------//

void Comm::finalize()
{
  MPI_Finalize();
}

//---------------------------------------------------------------------------//
// QUERY FUNCTIONS
//---------------------------------------------------------------------------//

int Comm::rank()
{
  int tmp = 0;
  MPI_Comm_rank(communicator, &tmp);

  // Postconditions
  Ensure(tmp >= 0);
  return tmp;
}

int Comm::size()
{
  int tmp = 0;
  MPI_Comm_size(communicator, &tmp);

  // Postconditions
  Ensure(tmp > 0);
  return tmp;
}

} // end namespace serment_comm

//---------------------------------------------------------------------------//
//              end of file MPI.cc
//---------------------------------------------------------------------------//
