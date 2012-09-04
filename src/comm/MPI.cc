//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MPI.cc
 * \brief  MPI 
 * \author Jeremy Roberts
 * \date   Aug 21, 2012
 */
//---------------------------------------------------------------------------//

#include "Comm.hh"

#include <vector>

#ifdef SERMENT_ENABLE_MPI

namespace serment_comm
{

// Set the default communicators
Communicator_t world          = MPI_COMM_WORLD;
Communicator_t global         = world;
Communicator_t local          = world;
Communicator_t communicator   = world;

// Initialize private static variables
bool   Comm::d_is_comm_built = false;
bool   Comm::d_is_global     = true;
double Comm::d_time = 0.0;
int    Comm::g_rank = 0;
//---------------------------------------------------------------------------//
// SETUP FUNCTIONS
//---------------------------------------------------------------------------//

void Comm::initialize(int &argc, char **&argv)
{
  int result = MPI_Init(&argc, &argv);
  d_time = 0.0;
  g_rank = Comm::rank();
  Ensure(result == MPI_SUCCESS);
}

void Comm::finalize()
{
  MPI_Finalize();
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
  Insist(communicator == world,
    "The current communicator must be world when setting up comm.");

  // Number of processes in world
  int world_size = Comm::size();

  Insist(N <= world_size,
    "The number of local groups cannot exceed the number of world processes");

  Insist(N > 0,
    "The number of local groups must be > 0");

  // Number of process per local communicator and remainder.
  int local_size = world_size / N;
  int remainder  = world_size - N * local_size;

  // Assign the number per local communicator.  Extra processes are put
  // in the first remainder processes.
  std::vector<int> proc_per_local(N, local_size);
  for (int i = 0; i < remainder; i++)
    proc_per_local[i] += 1;

  // Initialize the exclude list for the global communicator.  Processes
  // that are roots of local communicators will *not* be excluded.
  std::vector<int> exclude(world_size - N, 0);

  // Define the local color.
  int color = 0;
  int proc = 0;
  int exclude_index = 0;
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < proc_per_local[i]; j++, proc++)
    {
      if (proc == Comm::rank())
      {
        // This process gets the appropriate group color
        color = i;
      }
      if (j > 0)
      {
        // Exclude a process from the global communicator
        exclude[exclude_index++] = proc;
        if (proc == Comm::rank()) d_is_global = false;
      }
    }
  }

  int ierr;

  // Create the local communicator
  ierr = MPI_Comm_split(world, color, color, &local);

  // Create the global communicator
  MPI_Group world_g, global_g;
  ierr = MPI_Comm_group(world, &world_g);
  ierr = MPI_Group_excl(world_g, world_size - N, &exclude[0], &global_g);
  ierr = MPI_Comm_create(world, global_g, &global);
  ierr = MPI_Group_free(&global_g);
  ierr = MPI_Group_free(&world_g);

  d_is_comm_built = true;
  Ensure(!ierr);
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

#endif

//---------------------------------------------------------------------------//
//              end of file MPI.cc
//---------------------------------------------------------------------------//
