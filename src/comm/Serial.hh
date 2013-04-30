//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Serial.hh
 *  @brief  Serial comm implementation inline member definitions.
 *  @author Jeremy Roberts
 *  @date   Aug 21, 2012
 */
//---------------------------------------------------------------------------//

#ifndef serment_comm_SERIAL_HH_
#define serment_comm_SERIAL_HH_

#include "utilities/DBC.hh"
#include <ctime>

namespace serment_comm
{

//---------------------------------------------------------------------------//
// MPI Communicator
//---------------------------------------------------------------------------//

typedef int Communicator_t;

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
void Comm::set(const C &new_communicator)
{
  /* ... */
}

/// Free a communicator
inline void Comm::free()
{
  /* ... */
}

//---------------------------------------------------------------------------//
// BLOCKING SEND/RECEIVE OPERATIONS
//---------------------------------------------------------------------------//

template<class T>
int Comm::send(const T *buffer,
               int      size,
               int      destination,
               int      tag)
{
  return COMM_SUCCESS;
}

template<class T>
int Comm::receive(T   *buffer,
                  int  size,
                  int  source,
                  int  tag)
{
  return COMM_SUCCESS;
}

//---------------------------------------------------------------------------//
// BROADCAST
//---------------------------------------------------------------------------//

template<class T>
int Comm::broadcast(T  *buffer,
                    int size,
                    int root)
{
  return COMM_SUCCESS;
}

template<class T>
int Comm::broadcast(detran_utilities::SP<T> &buffer,
                    int                      root)
{
  return COMM_SUCCESS;
}

//---------------------------------------------------------------------------//
// REDUCTIONS
//---------------------------------------------------------------------------//

template<class T>
inline void Comm::sum(T &x, int to_node)
{
  /* ... */
}

template<class T>
inline void Comm::prod(T &x, int to_node)
{
  /* ... */
}

template<class T>
inline void Comm::min(T &x, int to_node)
{
  /* ... */
}

template<class T>
inline void Comm::max(T &x, int to_node)
{
  /* ... */
}

template<class T>
inline void Comm::sum(T *x, int n, int to_node)
{
  Require(x);
}

template<class T>
inline void Comm::prod(T  *x, int n, int to_node)
{
  Require(x);
}

template<class T>
inline void Comm::min(T *x, int n, int to_node)
{
  Require(x);
}

template<class T>
inline void Comm::max(T *x, int n, int to_node)
{
  Require(x);
}

//---------------------------------------------------------------------------//
// GLOBAL REDUCTIONS
//---------------------------------------------------------------------------//

template<class T>
inline void Comm::global_sum(T &x)
{
  /* ... */
}

template<class T>
inline void Comm::global_prod(T &x)
{
  /* ... */
}

template<class T>
inline void Comm::global_min(T &x)
{
  /* ... */
}

template<class T>
inline void Comm::global_max(T &x)
{
  /* ... */
}

template<class T>
inline void Comm::global_sum(T *x, int n)
{
  Require(x);
}

template<class T>
inline void Comm::global_prod(T  *x, int n)
{
  Require(x);
}

template<class T>
inline void Comm::global_min(T *x, int n)
{
  Require(x);
}

template<class T>
inline void Comm::global_max(T *x, int n)
{
  Require(x);
}

//---------------------------------------------------------------------------//
// BARRIER
//---------------------------------------------------------------------------//

inline void Comm::global_barrier()
{
  /* ... */
}

//---------------------------------------------------------------------------//
// TIMING \todo Need serial timer
//---------------------------------------------------------------------------//

inline void Comm::tic()
{
  d_time = (double) std::clock() / (double)CLOCKS_PER_SEC;
}

inline double Comm::toc()
{
  return (double) std::clock() / (double)CLOCKS_PER_SEC - d_time;
}

//---------------------------------------------------------------------------//
// UTILITIES
//---------------------------------------------------------------------------//

inline void Comm::partition(unsigned int &global_count,
                            unsigned int &local_start,
                            unsigned int &local_count)
{
  local_start = 0;
  local_count = global_count;
}

} // end namespace serment_comm

#endif // serment_comm_SERIAL_HH_

//---------------------------------------------------------------------------//
//              end of file Serial.hh
//---------------------------------------------------------------------------//
