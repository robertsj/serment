//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   LeakageOperator.hh
 *  @brief  LeakageOperator
 *  @author Jeremy Roberts
 *  @date   Aug 24, 2012
 */
//---------------------------------------------------------------------------//

#ifndef erme_LEAKAGEOPERATOR_HH_
#define erme_LEAKAGEOPERATOR_HH_

#include "ResponseOperator.hh"
#include "linear_algebra/Matrix.hh"

namespace erme
{

/**
 *  @class LeakageOperator
 *  @brief Leakage operator
 */
/**
 *  @example erme/test/test_LeakageOperator
 *
 *  Test of LeakageOperator class
 */
class LeakageOperator: public ResponseOperator,
                       public linear_algebra::Matrix
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<LeakageOperator>     SP_leakage;
  typedef linear_algebra::Vector::SP_vector         SP_vector;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param nodes    Pointer to node list
   *  @param indexer  Pointer to response indexer
   *  @param server   Pointer to response server
   */
  LeakageOperator(SP_nodelist nodes, SP_indexer indexer, SP_server server);

  /// Update the response matrix data
  void update();

  /**
   *  @brief Compute the net global leakage given a global moments vector
   *  @param x    Global moments vector
   */
  double leakage(linear_algebra::Vector &x);

  /// Display the global leakage vector
  void display_leakage();

private:

  //-------------------------------------------------------------------------//
  // PRIVATE DATA
  //-------------------------------------------------------------------------//

  /// Dotted with L*J, this gives the net leakage at global boundaries
  linear_algebra::Vector d_global_leakage;

  /// Vector for storing L*J
  linear_algebra::Vector d_L_times_moments;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

};

} // end namespace erme


#endif // erme_LEAKAGEOPERATOR_HH_

//---------------------------------------------------------------------------//
//              end of file LeakageOperator.hh
//---------------------------------------------------------------------------//
