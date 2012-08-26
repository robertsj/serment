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

namespace serment_comm
{

//---------------------------------------------------------------------------//
// MPI Communicator
//---------------------------------------------------------------------------//

extern MPI_Comm default_communicator;
extern MPI_Comm communicator;

typedef MPI_Comm Communicator_t;

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
// GLOBAL REDUCTIONS
//---------------------------------------------------------------------------//

template<class T>
inline void global_sum(T &x)
{
  // Copy data into send buffer
  T y = x;
  // Do global MPI reduction (result is on all processors) into x
  MPI_Allreduce(&y, &x, 1, MPI_Traits<T>::element_type(),
                MPI_SUM, communicator);
}

template<class T>
inline void global_prod(T &x)
{
  // copy data into send buffer
  T y = x;
  // Do global MPI reduction (result is on all processors) into x
  MPI_Allreduce(&y, &x, 1, MPI_Traits<T>::element_type(),
                MPI_PROD, communicator);
}

template<class T>
inline void global_min(T &x)
{
  // copy data into send buffer
  T y = x;
  // Do global MPI reduction (result is on all processors) into x
  MPI_Allreduce(&y, &x, 1, MPI_Traits<T>::element_type(),
                MPI_MIN, communicator);
}

template<class T>
inline void global_max(T &x)
{
  // copy data into send buffer
  T y = x;
  // Do global MPI reduction (result is on all processors) into x
  MPI_Allreduce(&y, &x, 1, MPI_Traits<T>::element_type(),
                MPI_MAX, communicator);
}

template<class T>
inline void global_sum(T *x, int n)
{
  Require (x);
  // Copy data into a send buffer
  std::vector<T> send_buffer(x, x + n);
  // Element-wise global reduction (result is on all processors) into x
  MPI_Allreduce(&send_buffer[0], x, n, MPI_Traits<T>::element_type(),
                MPI_SUM, communicator);
}

template<class T>
inline void global_prod(T  *x, int n)
{
  Require (x);
  // Copy data into a send buffer
  std::vector<T> send_buffer(x, x + n);
  // Element-wise global reduction (result is on all processors) into x
  MPI_Allreduce(&send_buffer[0], x, n, MPI_Traits<T>::element_type(),
                MPI_PROD, communicator);
}

template<class T>
inline void global_min(T *x, int n)
{
  Require (x);
  // Copy data into a send buffer
  std::vector<T> send_buffer(x, x + n);
  // Element-wise global reduction (result is on all processors) into x
  MPI_Allreduce(&send_buffer[0], x, n, MPI_Traits<T>::element_type(),
                MPI_MIN, communicator);
}

template<class T>
void global_max(T *x, int n)
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
