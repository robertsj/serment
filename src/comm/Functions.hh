//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Functions.hh
 * \brief  Defines the parallel communications interface.
 * \author Jeremy Roberts
 * \date   Aug 21, 2012
 * \node   Largely unchanged from Denovo source.
 */
//---------------------------------------------------------------------------//

#ifndef FUNCTIONS_HH_
#define FUNCTIONS_HH_

// Serment Comm
#include "Comm_Traits.hh"
#include "Global.hh"

namespace serment_comm
{

//---------------------------------------------------------------------------//
// SETUP FUNCTIONS
//---------------------------------------------------------------------------//

/*!
 * \brief Initialize a parallel job.
 */
void initialize(int &argc, char **&argv);

//---------------------------------------------------------------------------//
/*!
 * \brief Finish a parallel job.
 */
void finalize();


//---------------------------------------------------------------------------//
// QUERY FUNCTIONS
//---------------------------------------------------------------------------//

/*!
 * \brief Get the rank of the current processor.
 *
 * The rank is determined by the current communicator.
 */
int rank();

/*!
 * \brief Get the number of processors used for this job.
 *
 * The number of processes is determined by the current communicator.
 */
int size();

//---------------------------------------------------------------------------//
// BLOCKING SEND/RECEIVE OPERATIONS
//---------------------------------------------------------------------------//

/*!
 * \brief Do a point-to-point, blocking send.
 * \param buffer
 * \param size
 * \param destination
 * \param tag
 */
template<class T>
int send(const T *buffer, int size, int destination,
         int tag = Comm_Traits<T*>::tag);

/*!
 * \brief Do a point-to-point, blocking receive.
 * \param buffer
 * \param size
 * \param source
 * \param tag
 */
template<class T>
int receive(T *buffer, int size, int source,
            int tag = Comm_Traits<T*>::tag);

//---------------------------------------------------------------------------//
// BROADCAST
//---------------------------------------------------------------------------//

/*!
 * \brief Do a root-to-all broadcast.
 * \param buffer
 * \param size
 * \param root
 */
template<class T>
int broadcast(T *buffer, int size, int root);

//---------------------------------------------------------------------------//
// BARRIER FUNCTIONS
//---------------------------------------------------------------------------//

/// Set a global barrier for the communicator.
void global_barrier();

} // end namespace serment_comm

#endif // FUNCTIONS_HH_ 

//---------------------------------------------------------------------------//
//              end of file Functions.hh
//---------------------------------------------------------------------------//
