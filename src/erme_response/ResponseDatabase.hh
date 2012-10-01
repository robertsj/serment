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
#include "ResponseSource.hh"
#include "erme_geometry/DummyNode.hh"
#include "hdf5.h"

namespace erme_response
{

/**
 *  @class ResponseDatabase
 *  @brief Provides precomputed responses stored on disk
 *
 *  T
 *
 */
class ResponseDatabase
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<ResponseDatabase>  SP_responsedatabase;
  typedef erme_geometry::Node::SP_node            SP_node;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param node   Pointer to a response database node
   */
  ResponseDatabase();

  ~ResponseDatabase();

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL SOURCES MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /// Compute a response for the requested incident index
  void compute(SP_response response, ResponseIndex index);

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// HDF5 file id
  hid_t d_file_id;

  /// HDF5 filename
  std::string d_filename;

  /// HDF5 is open
  bool d_open;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

};

} // end namespace erme_response

#endif // RESPONSEDATABASE_HH_ 

//---------------------------------------------------------------------------//
//              end of file ResponseDatabase.hh
//---------------------------------------------------------------------------//
