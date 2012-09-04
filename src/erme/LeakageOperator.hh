//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LeakageOperator.hh
 * \brief  LeakageOperator 
 * \author Jeremy Roberts
 * \date   Aug 24, 2012
 */
//---------------------------------------------------------------------------//

#ifndef LEAKAGEOPERATOR_HH_
#define LEAKAGEOPERATOR_HH_

#include "ResponseOperator.hh"
#include "linear_algebra/Matrix.hh"

#include "DBC.hh"
#include "SP.hh"

namespace erme
{

/*!
 *  \class LeakageOperator
 *  \brief Leakage operator
 */
/*!
 *  \example erme/test/test_LeakageOperator
 *
 *  Test of LeakageOperator classs
 */
class LeakageOperator: public linear_algebra::Matrix,
                       public ResponseOperator
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran::SP<LeakageOperator> SP_leakageoperator;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   *  \param nodes    Pointer to node list
   *  \param indexer  Pointer to response indexer
   *  \param server   Pointer to response server
   */
  LeakageOperator(SP_nodelist nodes, SP_indexer indexer, SP_server server);

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


#endif // LEAKAGEOPERATOR_HH_ 

//---------------------------------------------------------------------------//
//              end of file LeakageOperator.hh
//---------------------------------------------------------------------------//
