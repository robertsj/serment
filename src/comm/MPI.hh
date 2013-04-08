//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   MPI.hh
 *  @brief  MPI communications interface
 *  @author Jeremy Roberts
 *  @date   Aug 21, 2012
 */
//---------------------------------------------------------------------------//

#ifndef serment_comm_MPI_HH_
#define serment_comm_MPI_HH_

#include "MPI_Traits.hh"
#include "utilities/DBC.hh"
#include <mpi.h>
#include <vector>
#include <iostream>

namespace serment_comm
{

//---------------------------------------------------------------------------//
// MPI Communicator
//---------------------------------------------------------------------------//

typedef MPI_Comm Communicator_t;

/// All processes, i.e. MPI_COMM_WORLD
extern Communicator_t world;

/// Subset of world that includes processes solving global problem
extern Communicator_t global;

/**
 *  Partitioning of world such that each local root is part of
 *  global.  Processes in local participate as response sources, and
 *  the root process participates as a response server for the global
 *  solve.
 */
extern Communicator_t local;

/// Current communicator
extern Communicator_t communicator;

//---------------------------------------------------------------------------//
// COMMUNICATORS
//---------------------------------------------------------------------------//

/// Set the communicator
template<class C>
inline void Comm::set(const C &new_communicator)
{
  communicator = new_communicator;
}

/// Free the communicators
inline void Comm::free()
{
  communicator = world;
  int ierr = 0;
  if (d_is_comm_built)
  {
    ierr = MPI_Comm_free(&local);
    if (d_is_global)
      ierr = MPI_Comm_free(&global);
    d_is_comm_built = false;
  }
  Ensure(!ierr);
}

//---------------------------------------------------------------------------//
// BLOCKING SEND/RECEIVE OPERATIONS
//---------------------------------------------------------------------------//

template<class T>
inline int Comm::send(const T *buffer,
                      int      size,
                      int      destination,
                      int      tag)
{
  MPI_Send(const_cast<T *>(buffer), size, MPI_Traits<T>::element_type(),
           destination, tag, communicator);
  return COMM_SUCCESS;
}

template<class T>
inline int Comm::receive(T   *buffer,
            int  size,
            int  source,
            int  tag)
{
  int count = 0;

  // get a handle to the MPI_Status
  MPI_Status status;

  // do the blocking receive
  MPI_Recv(buffer, size, MPI_Traits<T>::element_type(), source, tag,
           communicator, &status);

  // get the count of received data
  MPI_Get_count(&status, MPI_Traits<T>::element_type(), &count);
  return count;
}

//---------------------------------------------------------------------------//
// BROADCAST
//---------------------------------------------------------------------------//

template<class T>
inline int Comm::broadcast(T  *buffer,
                           int size,
                           int root)
{
  int r = MPI_Bcast(buffer, size, MPI_Traits<T>::element_type(),
                    root, communicator);
  return r;
}

//---------------------------------------------------------------------------//
// REDUCTIONS
//---------------------------------------------------------------------------//

template<class T>
inline void Comm::sum(T &x, int to_node)
{
  // Copy data into send buffer
  T y = x;
  // Do global MPI reduction (result is on all processors) into x
  MPI_Reduce(&y, &x, 1, MPI_Traits<T>::element_type(),
             MPI_SUM, to_node, communicator);
}

template<class T>
inline void Comm::prod(T &x, int to_node)
{
  // copy data into send buffer
  T y = x;
  // Do global MPI reduction (result is on all processors) into x
  MPI_Reduce(&y, &x, 1, MPI_Traits<T>::element_type(),
             MPI_PROD, to_node, communicator);
}

template<class T>
inline void Comm::min(T &x, int to_node)
{
  // copy data into send buffer
  T y = x;
  // Do global MPI reduction (result is on all processors) into x
  MPI_Reduce(&y, &x, 1, MPI_Traits<T>::element_type(),
             MPI_MIN, to_node, communicator);
}

template<class T>
inline void Comm::max(T &x, int to_node)
{
  // copy data into send buffer
  T y = x;
  // Do global MPI reduction (result is on all processors) into x
  MPI_Reduce(&y, &x, 1, MPI_Traits<T>::element_type(),
             MPI_MAX, to_node, communicator);
}

template<class T>
inline void Comm::sum(T *x, int n, int to_node)
{
  Require (x);
  // Copy data into a send buffer
  std::vector<T> send_buffer(x, x + n);
  // Element-wise global reduction (result is on all processors) into x
  MPI_Reduce(&send_buffer[0], x, n, MPI_Traits<T>::element_type(),
             MPI_SUM, to_node, communicator);
}

template<class T>
inline void Comm::prod(T  *x, int n, int to_node)
{
  Require (x);
  // Copy data into a send buffer
  std::vector<T> send_buffer(x, x + n);
  // Element-wise global reduction (result is on all processors) into x
  MPI_Reduce(&send_buffer[0], x, n, MPI_Traits<T>::element_type(),
             MPI_PROD, to_node, communicator);
}

template<class T>
inline void Comm::min(T *x, int n, int to_node)
{
  Require (x);
  // Copy data into a send buffer
  std::vector<T> send_buffer(x, x + n);
  // Element-wise global reduction (result is on all processors) into x
  MPI_Reduce(&send_buffer[0], x, n, MPI_Traits<T>::element_type(),
             MPI_MIN, to_node, communicator);
}

template<class T>
inline void Comm::max(T *x, int n, int to_node)
{
  Require (x);
  // Copy data into a send buffer
  std::vector<T> send_buffer(x, x + n);
  // Element-wise global reduction (result is on all processors) into x
  MPI_Reduce(&send_buffer[0], x, n, MPI_Traits<T>::element_type(),
             MPI_MAX, to_node, communicator);
}

//---------------------------------------------------------------------------//
// GLOBAL REDUCTIONS
//---------------------------------------------------------------------------//

template<class T>
inline void Comm::global_sum(T &x)
{
  // Copy data into send buffer
  T y = x;
  // Do global MPI reduction (result is on all processors) into x
  MPI_Allreduce(&y, &x, 1, MPI_Traits<T>::element_type(),
                MPI_SUM, communicator);
}

template<class T>
inline void Comm::global_prod(T &x)
{
  // copy data into send buffer
  T y = x;
  // Do global MPI reduction (result is on all processors) into x
  MPI_Allreduce(&y, &x, 1, MPI_Traits<T>::element_type(),
                MPI_PROD, communicator);
}

template<class T>
inline void Comm::global_min(T &x)
{
  // copy data into send buffer
  T y = x;
  // Do global MPI reduction (result is on all processors) into x
  MPI_Allreduce(&y, &x, 1, MPI_Traits<T>::element_type(),
                MPI_MIN, communicator);
}

template<class T>
inline void Comm::global_max(T &x)
{
  // copy data into send buffer
  T y = x;
  // Do global MPI reduction (result is on all processors) into x
  MPI_Allreduce(&y, &x, 1, MPI_Traits<T>::element_type(),
                MPI_MAX, communicator);
}

template<class T>
inline void Comm::global_sum(T *x, int n)
{
  Require (x);
  // Copy data into a send buffer
  std::vector<T> send_buffer(x, x + n);
  // Element-wise global reduction (result is on all processors) into x
  MPI_Allreduce(&send_buffer[0], x, n, MPI_Traits<T>::element_type(),
                MPI_SUM, communicator);
}

template<class T>
inline void Comm::global_prod(T  *x, int n)
{
  Require (x);
  // Copy data into a send buffer
  std::vector<T> send_buffer(x, x + n);
  // Element-wise global reduction (result is on all processors) into x
  MPI_Allreduce(&send_buffer[0], x, n, MPI_Traits<T>::element_type(),
                MPI_PROD, communicator);
}

template<class T>
inline void Comm::global_min(T *x, int n)
{
  Require (x);
  // Copy data into a send buffer
  std::vector<T> send_buffer(x, x + n);
  // Element-wise global reduction (result is on all processors) into x
  MPI_Allreduce(&send_buffer[0], x, n, MPI_Traits<T>::element_type(),
                MPI_MIN, communicator);
}

template<class T>
inline void Comm::global_max(T *x, int n)
{
  Require (x);
  // Copy data into a send buffer
  std::vector<T> send_buffer(x, x + n);
  // Element-wise global reduction (result is on all processors) into x
  MPI_Allreduce(&send_buffer[0], x, n, MPI_Traits<T>::element_type(),
                MPI_MAX, communicator);
}

//---------------------------------------------------------------------------//
// BARRIER
//---------------------------------------------------------------------------//

inline void Comm::global_barrier()
{
  MPI_Barrier(communicator);
}

//---------------------------------------------------------------------------//
// TIMING
//---------------------------------------------------------------------------//

inline void Comm::tic()
{
  d_time = MPI_Wtime();
}

inline double Comm::toc()
{
  return MPI_Wtime() - d_time;
}

//---------------------------------------------------------------------------//
// UTILITIES
//---------------------------------------------------------------------------//

/*
 *  This implementation assigns each process an equal size chunk of the
 *  total count.  Extra counts are assigned one-per-process starting with
 *  the last process in the communicator.
 *
 */
inline void Comm::partition(unsigned int &global_count,
                            unsigned int &local_start,
                            unsigned int &local_count)
{

  // Number of things per process
  std::vector<unsigned int> number_per_process(Comm::size(), 0);

  // Master
  if (Comm::rank() == 0)
  {
    // Preconditions
    Insist(global_count > 0,
      "Global count must be a positive in order to partition it");

    // Initial guess for count per process and the remainder.
    int number_per_process_guess = global_count / Comm::size();
    int remainder = global_count - number_per_process_guess * Comm::size();

    // Assign the number per process, putting the extras on
    // the processes in reverse.  If there are more processes
    // than nodes, put the things on the first processes.
    int number_processes = Comm::size();
    if (!number_per_process_guess)
    {
      number_processes = global_count;
      number_per_process_guess = 1;
      remainder = 0;
    }

    for (int i = 0; i < number_processes; i++)
      number_per_process[i] = number_per_process_guess;

    for (int i = 1; i <= remainder; i++)
      number_per_process[Comm::size() - i] += 1;
  }

  // Broadcast the total number of things
  Comm::broadcast(&global_count, 1, 0);

  // Broadcast the number of nodes per process.
  Comm::broadcast(&number_per_process[0], number_per_process.size(), 0);

  // Local start and size
  local_start = 0;
  for (int i = 0; i < Comm::rank(); i++)
    local_start += number_per_process[i];
  local_count = number_per_process[Comm::rank()];

}

} // end namespace serment_comm

#endif // serment_comm_MPI_HH_

//---------------------------------------------------------------------------//
//              end of file MPI.hh
//---------------------------------------------------------------------------//
