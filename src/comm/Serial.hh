//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Serial.hh
 * \brief  Serial comm implementation inline member definitions.
 * \author Jeremy Roberts
 * \date   Aug 21, 2012
 */
//---------------------------------------------------------------------------//

#ifndef SERIAL_HH_
#define SERIAL_HH_

// Detran utilities
#include "DBC.hh"

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
template<class C>
void Comm::free(C &new_communicator)
{
  /* ... */
}

//---------------------------------------------------------------------------//
// BLOCKING SEND/RECEIVE OPERATIONS
//---------------------------------------------------------------------------//

template<class T>
int send(const T *buffer,
         int      size,
         int      destination,
         int      tag)
{
  return COMM_SUCCESS;
}

template<class T>
int receive(T   *buffer,
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
int broadcast(T  *buffer,
              int size,
              int root)
{
  return COMM_SUCCESS;
}

//---------------------------------------------------------------------------//
// REDUCTIONS
//---------------------------------------------------------------------------//

template<class T>
inline void sum(T &x, int to_node)
{
  return x;
}

template<class T>
inline void prod(T &x, int to_node)
{
  return x;
}

template<class T>
inline void min(T &x, int to_node)
{
  /* ... */
}

template<class T>
inline void max(T &x, int to_node)
{
  /* ... */
}

template<class T>
inline void sum(T *x, int n, int to_node)
{
  Require(x);
}

template<class T>
inline void prod(T  *x, int n, int to_node)
{
  Require(x);
}

template<class T>
inline void min(T *x, int n, int to_node)
{
  Require(x);
}

template<class T>
inline void max(T *x, int n, int to_node)
{
  Require(x);
}

//---------------------------------------------------------------------------//
// GLOBAL REDUCTIONS
//---------------------------------------------------------------------------//

template<class T>
inline void global_sum(T &x)
{
  return x;
}

template<class T>
inline void global_prod(T &x)
{
  return x;
}

template<class T>
inline void global_min(T &x)
{
  /* ... */
}

template<class T>
inline void global_max(T &x)
{
  /* ... */
}

template<class T>
inline void global_sum(T *x, int n)
{
  Require(x);
}

template<class T>
inline void global_prod(T  *x, int n)
{
  Require(x);
}

template<class T>
inline void global_min(T *x, int n)
{
  Require(x);
}

template<class T>
inline void global_max(T *x, int n)
{
  Require(x);
}

//---------------------------------------------------------------------------//
// BARRIER FUNCTIONS
//---------------------------------------------------------------------------//

void global_barrier()
{
  /* ... */
}

} // end namespace serment_comm

#endif // SERIAL_HH_ 

//---------------------------------------------------------------------------//
//              end of file Serial.hh
//---------------------------------------------------------------------------//
