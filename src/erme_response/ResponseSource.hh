//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseSource.hh
 * \brief  ResponseSource class definition
 * \author Jeremy Roberts
 * \date   Aug 28, 2012
 */
//---------------------------------------------------------------------------//

#ifndef RESPONSESOURCE_HH_
#define RESPONSESOURCE_HH_

#include "DBC.hh"
#include "SP.hh"
#include "NodeResponse.hh"

namespace erme_response
{

/*!
 *  \class ResponseSource
 *  \brief Abstract response source
 *
 *  A ResponseSource provides its ResponseServer with responses for use
 *  in the global solve.  Each ResponseSource is unique for a given Node.
 *  The ResponseSource represents the interface between Serment and local
 *  solvers such as Detran.  Each concrete Node implementation must have
 *  a corresponding concrete ResponseSource and
 *  ResponseSourceFactory::build specialization.
 *
 *
 */
class ResponseSource
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran::SP<ResponseSource>    SP_source;
  typedef unsigned int                  size_t;
  typedef NodeResponse::SP_response     SP_response;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   */
  ResponseSource(){};

  /// Virtual destructor
  virtual ~ResponseSource(){}

  void update(const double keff)
  {
    d_keff = keff;
  }

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE
  //-------------------------------------------------------------------------//

  /// Compute a response for the requested incident index
  virtual void compute(SP_response response, ResponseIndex index) = 0;

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// K-eigenvalue
  double d_keff;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//


};

} // end namespace erme_response

#endif // RESPONSESOURCE_HH_ 

//---------------------------------------------------------------------------//
//              end of file ResponseSource.hh
//---------------------------------------------------------------------------//
