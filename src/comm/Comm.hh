//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Comm.hh
 * \brief  Comm 
 * \author Jeremy Roberts
 * \date   Aug 21, 2012
 */
//---------------------------------------------------------------------------//

#ifndef COMM_HH_
#define COMM_HH_

#include "serment_config.h"

#include "Comm_Traits.hh"

namespace serment_comm
{

/*!
 *  \class Comm
 *  \brief Parallel communication interface
 *
 *  This is an easy API for using MPI (or some other parallel library).
 *
 */

class Comm
{

public:

  //---------------------------------------------------------------------------//
  // SETUP FUNCTIONS
  //---------------------------------------------------------------------------//

  /*!
   * \brief Initialize a parallel job.
   */
  static void initialize(int &argc, char **&argv);

  //---------------------------------------------------------------------------//
  /*!
   * \brief Finish a parallel job.
   */
  static void finalize();


  //---------------------------------------------------------------------------//
  // QUERY FUNCTIONS
  //---------------------------------------------------------------------------//

  /*!
   * \brief Get the rank of the current processor.
   *
   * The rank is determined by the current communicator.
   */
  static int rank();

  /*!
   * \brief Get the number of processors used for this job.
   *
   * The number of processes is determined by the current communicator.
   */
  static int size();

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
  static int send(const T *buffer, int size, int destination,
             int tag = Comm_Traits<T*>::tag);

  /*!
   * \brief Do a point-to-point, blocking receive.
   * \param buffer
   * \param size
   * \param source
   * \param tag
   */
  template<class T>
  static int receive(T *buffer, int size, int source,
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
  static int broadcast(T *buffer, int size, int root);

  //---------------------------------------------------------------------------//
  // BARRIER FUNCTIONS
  //---------------------------------------------------------------------------//

  /// Set a global barrier for the communicator.
  static void global_barrier();

  //---------------------------------------------------------------------------//
  // TIMING
  //---------------------------------------------------------------------------//

  /// Start the timer
  static void tic();

  /// Return time elapsed after tic() call.
  static double toc();

//private:

  /// Stored time information
  static double d_time;

};

} // end namespace serment_comm

#ifdef SERMENT_ENABLE_MPI
#include "MPI.hh"
#else
#include "Serial.hh"
#endif

#endif // COMM_HH_ 

//---------------------------------------------------------------------------//
//              end of file Comm.hh
//---------------------------------------------------------------------------//
