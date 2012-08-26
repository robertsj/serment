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

typedef int Communicator_t;

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
void global_max(T *x, int n)
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
