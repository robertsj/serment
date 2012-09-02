//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MPI.hh
 * \brief  MPI communications interface
 * \author Jeremy Roberts
 * \date   Aug 21, 2012
 */
//---------------------------------------------------------------------------//

#ifndef MPI_HH_
#define MPI_HH_

// Serment Comm
#include "MPI_Traits.hh"

// Detran utilities
#include "DBC.hh"

// System
#include <mpi.h>
#include <vector>
#include <iostream>

namespace serment_comm
{

//---------------------------------------------------------------------------//
// MPI Communicator
//---------------------------------------------------------------------------//

typedef MPI_Comm Communicator_t;

/*!
 *  All processes, i.e. MPI_COMM_WORLD
 */
extern Communicator_t world;

/*!
 *  Subset of world that includes processes solving global problem
 */
extern Communicator_t global;

/*!
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
  int ierr;
  if (d_is_comm_built)
  {
    ierr = MPI_Comm_free(&local);
    if (d_is_global)
      ierr = MPI_Comm_free(&global);
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

} // end namespace serment_comm

#endif // MPI_HH_ 

//---------------------------------------------------------------------------//
//              end of file MPI.hh
//---------------------------------------------------------------------------//
