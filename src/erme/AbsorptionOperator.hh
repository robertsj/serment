//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   AbsorptionOperator.hh
 * \brief  AbsorptionOperator
 * \author Jeremy Roberts
 * \date   Aug 24, 2012
 */
//---------------------------------------------------------------------------//

#ifndef ABSORPTIONOPERATOR_HH_
#define ABSORPTIONOPERATOR_HH_

#include "ResponseOperator.hh"
#include "linear_algebra/Vector.hh"

namespace erme
{

/*!
 *  \class AbsorptionOperator
 *  \brief Converts a global moments vector into a global absorption rate
 *
 */
/*!
 *  \example erme/test/test_AbsorptionOperator
 *
 *  Test of AbsorptionOperator class
 */
class AbsorptionOperator: public linear_algebra::Vector,
                          public ResponseOperator
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran::SP<AbsorptionOperator>       SP_fission;


  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   *  \param nodes    Pointer to node list
   *  \param indexer  Pointer to response indexer
   *  \param server   Pointer to response server
   */
  AbsorptionOperator(SP_nodelist nodes, SP_indexer indexer, SP_server server);

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

#endif // ABSORPTIONOPERATOR_HH_ 

//---------------------------------------------------------------------------//
//              end of file AbsorptionOperator.hh
//---------------------------------------------------------------------------//
