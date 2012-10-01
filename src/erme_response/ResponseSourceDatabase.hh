//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ResponseSourceDatabase.hh
 *  @brief  ResponseSourceDatabase
 *  @author Jeremy Roberts
 *  @date   Sep 30, 2012
 */
//---------------------------------------------------------------------------//

#ifndef erme_response_RESPONSESOURCEDATABASE_HH_
#define erme_response_RESPONSESOURCEDATABASE_HH_

#include "ResponseSource.hh"
#include "ResponseDatabase.hh"

namespace erme_response
{

/**
 *  @class ResponseDatabase
 *  @brief Provides precomputed responses from a database
 */
class ResponseSourceDatabase: public ResponseSource
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef ResponseDatabase::SP_rfdb             SP_rfdb;
  typedef erme_geometry::Node::SP_node          SP_node;

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

  /// Response function database
  SP_rfdb d_rfdb;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

};

} // end namespace erme_response

#endif // erme_response_RESPONSESOURCEDATABASE_HH_

//---------------------------------------------------------------------------//
//              end of file ResponseSourceDatabase.hh
//---------------------------------------------------------------------------//
