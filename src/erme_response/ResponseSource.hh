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

namespace erme_response
{

/*!
 *  \class ResponseSource
 *  \brief Abstract response source
 *
 *  The ResponseServer supplies nodal responses to various
 *  clients, and a concrete implementation of a ResponseSource
 *  provides the responses to the server.  The source can be
 *  an on-the-fly generation of responses or database.
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

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   */
  ResponseSource();

  /// Virtual destructor
  virtual ~ResponseSource(){}

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE
  //-------------------------------------------------------------------------//




private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//


};

} // end namespace erme_response

#endif // RESPONSESOURCE_HH_ 

//---------------------------------------------------------------------------//
//              end of file ResponseSource.hh
//---------------------------------------------------------------------------//
