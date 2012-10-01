//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseSourceDatabase.hh
 * \brief  ResponseSourceDatabase 
 * \author Jeremy Roberts
 * \date   Sep 30, 2012
 */
//---------------------------------------------------------------------------//

#ifndef erme_response_RESPONSESOURCEDATABASE_HH_
#define erme_response_RESPONSESOURCEDATABASE_HH_

namespace erme_response
{

/**
 *  @class ResponseDatabase
 *  @brief Provides precomputed responses stored on disk
 *
 *  T
 *
 */
class ResponseSourceDatabase: public ResponseSource
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<ResponseSourceDatabase>  SP_responsedatabase;
  typedef erme_geometry::Node::SP_node                  SP_node;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param node   Pointer to a response database node
   */
  ResponseSourceDatabase(SP_node node);

  /// Virtual destructor
  virtual ~ResponseSourceDatabase();

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

#endif // erme_response_RESPONSESOURCEDATABASE_HH_

//---------------------------------------------------------------------------//
//              end of file ResponseSourceDatabase.hh
//---------------------------------------------------------------------------//
