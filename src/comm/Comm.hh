//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Comm.hh
 *  @brief  Comm
 *  @author Jeremy Roberts
 *  @date   Aug 21, 2012
 */
//---------------------------------------------------------------------------//

#ifndef serment_comm_COMM_HH_
#define serment_comm_COMM_HH_

// Configuration
#include "serment_config.h"

// Serment Comm
#include "Comm_Traits.hh"

namespace serment_comm
{

/**
 *  @page Serment Parallel Structure
 *
 *  Serment employs a flexible parallel decomposition that should adapt
 *  well to just about any architecture.  The decomposition largely
 *  follows the physical decomposition implied by the method in which
 *  a global domain (e.g. a reactor core) is separated into distinct
 *  nodes (e.g. assemblies) that are linked via approximate boundary
 *  conditions.
 *
 *  Suppose we have a compute cluster with \f$ N \f$ processes
 *  available.  The set of all these processes constitutes the
 *  @e world communicator.
 *
 *  Since the bulk of the computation occurs during computation of
 *  the nodal boundary conditions, the number of processes that need
 *  to participate in the global solution might be a small subset of
 *  the total number of process available.  Consequently, we split
 *  @e world into \f$ M \f$ new communicators denoted \e local.
 *
 *  The root process of each \e local communicator is then included
 *  in a \e global communicator.  Processes in \e global participate
 *  in the global problem solution, i.e. they solve the response
 *  matrix equations.
 *
 *  Within a \e local communicator, the root process acts as the
 *  response server (to the global problem), and all processes
 *  serve as response sources (i.e. they generate response functions
 *  as needed by the server).
 *
 *  It should be noted an implicit third level of decomposition
 *  exists on systems with many core architectures.  Each process
 *  in a \e local communicator might reside on a physical compute
 *  node with several CPU cores or a GP-GPU accelerator.  In this
 *  case, the computation of a single response can also be
 *  parallelized via threading.
 *
 *  This scheme is flexible, and so it should allow a user to tailor
 *  a decomposition for maximizing a particular system's capabilities.
 *
 */


/**
 *  @class Comm
 *  @brief Parallel communication interface
 *
 *  This is an easy API for using MPI (or some other parallel library).
 *
 *  @todo Consider implementing Comm as a singleton so that better
 *        encapsulation can be used
 */
class Comm
{

public:

  //---------------------------------------------------------------------------//
  // SETUP FUNCTIONS
  //---------------------------------------------------------------------------//

  /// Initialize a parallel job.
  static void initialize(int &argc, char **&argv);

  /// Finish a parallel job.
  static void finalize();

  //---------------------------------------------------------------------------//
  // COMMUNICATORS
  //---------------------------------------------------------------------------//

  /**
   *  @brief Create communicators
   *  @param N  Number of local groups to create
   */
  static void setup_communicators(const unsigned int N = 1);

  /// Set the communicator
  template<class C>
  static void set(const C &new_communicator);

  /// Free a communicator
  static void free();

  /// Am I part of the global communicator?
  static bool is_global()
  {
    return d_is_global;
  }

  /// Have comm groups been built?  This could be put
  /// into the init routine, but then a new parameter would be added.
  static bool is_comm_built()
  {
    return d_is_comm_built;
  }

  //---------------------------------------------------------------------------//
  // QUERY FUNCTIONS
  //---------------------------------------------------------------------------//

  /**
   *  @brief Get the rank of the current processor.
   *
   * The rank is determined by the current communicator.
   */
  static int rank();

  /**
   *  @brief Get the number of processors used for this job.
   *
   * The number of processes is determined by the current communicator.
   */
  static int size();

  /// Return the rank of a process within the world communicator
  static int world_rank()
  {
    return d_world_rank;
  }

  /// Return the local group of a process
  static int local_group()
  {
    return d_local_group;
  }

  //---------------------------------------------------------------------------//
  // BLOCKING SEND/RECEIVE OPERATIONS
  //---------------------------------------------------------------------------//

  /**
   *  @brief Do a point-to-point, blocking send.
   *  @param buffer
   *  @param size
   *  @param destination
   *  @param tag
   */
  template<class T>
  static int send(const T *buffer, int size, int destination,
             int tag = Comm_Traits<T*>::tag);

  /**
   *  @brief Do a point-to-point, blocking receive.
   *  @param buffer
   *  @param size
   *  @param source
   *  @param tag
   */
  template<class T>
  static int receive(T *buffer, int size, int source,
                     int tag = Comm_Traits<T*>::tag);

  //---------------------------------------------------------------------------//
  // BROADCAST
  //---------------------------------------------------------------------------//

  /**
   *  @brief Do a root-to-all broadcast.
   *  @param buffer
   *  @param size
   *  @param root
   */
  template<class T>
  static int broadcast(T *buffer, int size, int root = 0);

  //---------------------------------------------------------------------------//
  // REDUCTIONS
  //---------------------------------------------------------------------------//

  /// Sum of a scalar variable.
  template <class T>
  static void sum(T &x, int to_node);

  /// Global product of a scalar variable.
  template <class T>
  static void prod(T &x, int to_node);

  /// Global minimum of a scalar variable.
  template <class T>
  static void min(T &x, int to_node);

  /// Global maximum of a scalar variable.
  template <class T>
  static void max(T &x, int to_node);

  /// Element-wise, global sum of an array.
  template <class T>
  static void sum(T *x, int n, int to_node);

  /// Element-wise, global product of an array.
  template <class T>
  static void prod(T *x, int n, int to_node);

  /// Element-wise, global minimum of an array.
  template <class T>
  static void min(T *x, int n, int to_node);

  /// Element-wise, global maximum of an array.
  template <class T>
  static void max(T *x, int n, int to_node);

  //---------------------------------------------------------------------------//
  // GLOBAL REDUCTIONS
  //---------------------------------------------------------------------------//

  /// Global sum of a scalar variable.
  template <class T>
  static void global_sum(T &x);

  /// Global product of a scalar variable.
  template <class T>
  static void global_prod(T &x);

  /// Global minimum of a scalar variable.
  template <class T>
  static void global_min(T &x);

  /// Global maximum of a scalar variable.
  template <class T>
  static void global_max(T &x);

  /// Element-wise, global sum of an array.
  template <class T>
  static void global_sum(T *x, int n);

  /// Element-wise, global product of an array.
  template <class T>
  static void global_prod(T *x, int n);

  /// Element-wise, global minimum of an array.
  template <class T>
  static void global_min(T *x, int n);

  /// Element-wise, global maximum of an array.
  template <class T>
  static void global_max(T *x, int n);

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

  //---------------------------------------------------------------------------//
  // UTILITY
  //---------------------------------------------------------------------------//

  /**
   *  @brief Partition a number represented an amount of work, etc.
   *  @param global_count   Total number of things to be partitioned
   *  @param local_start    My local index in array of things
   *  @param local_count    Number of local things I own
   */
  static void partition(unsigned int &global_count,
                        unsigned int &local_start,
                        unsigned int &local_count);

private:

  /// Are the communicators built?
  static bool d_is_comm_built;
  /// Am I on the global comm?
  static bool d_is_global;
  /// Stored time information
  static double d_time;
  /// Rank in world
  static int d_world_rank;
  /// Local group (= to color)
  static int d_local_group;

};

} // end namespace serment_comm

#ifdef SERMENT_ENABLE_MPI
#include "MPI.hh"
#else
#include "Serial.hh"
#endif

#endif // serment_comm_COMM_HH_

//---------------------------------------------------------------------------//
//              end of file Comm.hh
//---------------------------------------------------------------------------//
