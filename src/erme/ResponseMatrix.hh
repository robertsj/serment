//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ResponseMatrix.hh
 *  @brief  ResponseMatrix class definition
 *  @author Jeremy Roberts
 *  @date   Aug 24, 2012
 */
//---------------------------------------------------------------------------//

#ifndef erme_RESPONSEMATRIX_HH_
#define erme_RESPONSEMATRIX_HH_

#include "ResponseOperator.hh"
#include "linear_algebra/Matrix.hh"

namespace erme
{

/**
 *  @class ResponseMatrix
 *  @brief Response matrix operator
 */
/**
 *  @example erme/test/test_ResponseMatrix
 *
 *  Test of ResponseMatrix class
 */
class ResponseMatrix: public linear_algebra::Matrix,
                      public ResponseOperator
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<ResponseMatrix>      SP_responsematrix;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param nodes    Pointer to node list
   *  @param indexer  Pointer to response indexer
   *  @param server   Pointer to response server
   */
  ResponseMatrix(SP_nodelist nodes, SP_indexer indexer, SP_server server);

  /// Update the response matrix data
  void update();

private:

  //-------------------------------------------------------------------------//
  // PRIVATE DATA
  //-------------------------------------------------------------------------//


  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//



};

} // end namespace erme

#endif // erme_RESPONSEMATRIX_HH_

//---------------------------------------------------------------------------//
//              end of file ResponseMatrix.hh
//---------------------------------------------------------------------------//
