//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   FissionOperator.hh
 * \brief  FissionOperator 
 * \author Jeremy Roberts
 * \date   Aug 24, 2012
 */
//---------------------------------------------------------------------------//

#ifndef FISSIONOPERATOR_HH_
#define FISSIONOPERATOR_HH_

#include "ResponseOperator.hh"
#include "linear_algebra/Vector.hh"

namespace erme
{

/*!
 *  \class FissionOperator
 *  \brief Converts a global moments vector into a global fission rate
 *
 */
/*!
 *  \example erme/test/test_FissionOperator
 *
 *  Test of FissionOperator class
 */
class FissionOperator: public linear_algebra::Vector,
                       public ResponseOperator
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<FissionOperator>       SP_fission;


  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   *  \param nodes    Pointer to node list
   *  \param indexer  Pointer to response indexer
   *  \param server   Pointer to response server
   */
  FissionOperator(SP_nodelist nodes, SP_indexer indexer, SP_server server);

  /// Update the vector data
  void update();

private:

  //-------------------------------------------------------------------------//
  // PRIVATE DATA
  //-------------------------------------------------------------------------//


  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//



};

} // end namespace detran

#endif // FISSIONOPERATOR_HH_ 

//---------------------------------------------------------------------------//
//              end of file FissionOperator.hh
//---------------------------------------------------------------------------//
