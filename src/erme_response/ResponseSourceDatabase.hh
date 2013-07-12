//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ResponseSourceDatabase.hh
 *  @brief ResponseSourceDatabase
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

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

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef ResponseDatabase::SP_rfdb             SP_rfdb;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param node   Pointer to a response database node
   */
  ResponseSourceDatabase(SP_node node, SP_indexer indexer);

  /// Virtual destructor
  virtual ~ResponseSourceDatabase();

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL SOURCES MUST IMPLEMENT THESE
  //--------------------------------------------------------------------------//

  /// Compute a response for the requested incident index
  void compute(SP_response response, const ResponseIndex &index);

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Response function database
  SP_rfdb d_rfdb;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

};

} // end namespace erme_response

#endif // erme_response_RESPONSESOURCEDATABASE_HH_

//----------------------------------------------------------------------------//
//              end of file ResponseSourceDatabase.hh
//----------------------------------------------------------------------------//
