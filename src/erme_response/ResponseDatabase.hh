//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseDatabase.hh
 * \brief  ResponseDatabase 
 * \author Jeremy Roberts
 * \date   Sep 24, 2012
 */
//---------------------------------------------------------------------------//

#ifndef RESPONSEDATABASE_HH_
#define RESPONSEDATABASE_HH_

#include "NodeResponse.hh"
#include "ResponseIndex.hh"
#include "ioutils/IO_HDF5_Traits.hh"
#include <map>
#include <string>
#include <vector>
#include "hdf5.h"

namespace erme_response
{

/**
 *  @class ResponseDatabase
 *  @brief Provides precomputed responses stored on disk
 *
 *  This first implementation uses a singleton pattern.  The reason for
 *  this is so that just one instance of the database is created for
 *  all db-derived responses.  This model may break down when several
 *  local groups are used.
 *
 *  We take the following approach.  At the construction of the response
 *  server, the database is read in completely and then broadcasted to
 *  all nodes.  Interpolation and expansion happens using data in memory.
 *
 *  The routines for reading the data can alternatively be used to read
 *  data in as needed.  We'll leave that for later.
 */
class ResponseDatabase
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<ResponseDatabase>  SP_rfdb;
  typedef NodeResponse::SP_response               SP_response;
  typedef std::vector<SP_response>                vec_response;
  typedef struct
  {
    // number of keff terms (for interpolating/expanding)
    int number_keffs;
    // interpolate = 0, expand = 1
    int scheme;
    // vector of keffs
    std::vector<double> keffs;
    // keff in vacuum (for analytic summation of series)
    double kvac;
    // response data
    vec_response responses;
  } DBResponse;
  typedef std::map<std::string, DBResponse>       map_response;
  typedef map_response::iterator                  response_it;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief SP constructor
   *  @param filename  name of hdf5 file
   */
  static SP_rfdb Create(std::string filename);

  /// SP Accessor
  static SP_rfdb Instance();

  /// Destructor
  ~ResponseDatabase();

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /**
   *  @brief Get the response for a node, index, and keff
   *  @param nodename  name of node to look up in database
   *  @param response  response container to fill
   *  @param index     index of incident response
   *  @param keff      eigenvalue for requested response
   */
  void get(std::string     nodename,
           SP_response     response,
           ResponseIndex   index,
           const double    keff);


  std::string filename() const
  {
    return d_filename;
  }

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Single instance
  static SP_rfdb d_instance;

  /// HDF5 filename
  std::string d_filename;

  /// HDF5 is open
  bool d_open;

  /// HDF5 file id
  hid_t d_file_id;

  /// Number of nodes in the database
  hsize_t d_number_nodes;

  /**
   *  Response data
   *
   *  All response data is stored as a function of keff or as a function
   *  of expansion index for each node in the database.
   *
   *  ["nodename"][keff][data]
   */
  map_response d_responses;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param filename   name of hdf5 file
   */
  ResponseDatabase(std::string filename);

  /// Read a scalar (int or double) attribute
  template <class T>
  bool read_scalar_attribute(hid_t group, const char* name, T &value);


  /// Interpolate
  double interpolate(double x, double x0, double x1, double y0, double y1);

};

} // end namespace erme_response

#include "ResponseDatabase.i.hh"

#endif // RESPONSEDATABASE_HH_ 

//---------------------------------------------------------------------------//
//              end of file ResponseDatabase.hh
//---------------------------------------------------------------------------//
