//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseMatrix.hh
 * \brief  ResponseMatrix class definition
 * \author Jeremy Roberts
 * \date   Aug 24, 2012
 */
//---------------------------------------------------------------------------//

#ifndef RESPONSEMATRIX_HH_
#define RESPONSEMATRIX_HH_

#include "ResponseOperator.hh"
#include "linear_algebra/Matrix.hh"

#include "DBC.hh"
#include "SP.hh"

namespace erme
{

/*!
 *  \class ResponseMatrix
 *  \brief Response matrix operator
 */
class ResponseMatrix: public linear_algebra::Matrix,
                      public ResponseOperator
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran::SP<ResponseMatrix>                  SP_responsematrix;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   *  \param indexer  Pointer to response indexer
   *  \param server   Pointer to response server
   */
  ResponseMatrix(SP_nodelist nodes, SP_indexer indexer, SP_server server);

  /// Update the response matrix data for a new eigenvalue
  void update(const double keff);

private:

  //-------------------------------------------------------------------------//
  // PRIVATE DATA
  //-------------------------------------------------------------------------//


  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//



};

} // end namespace erme

#endif // RESPONSEMATRIX_HH_ 

//---------------------------------------------------------------------------//
//              end of file ResponseMatrix.hh
//---------------------------------------------------------------------------//
