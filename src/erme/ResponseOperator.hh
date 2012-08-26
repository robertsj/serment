//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseOperator.hh
 * \brief  ResponseOperator class definition
 * \author Jeremy Roberts
 * \date   Aug 24, 2012
 */
//---------------------------------------------------------------------------//

#ifndef RESPONSEOPERATOR_HH_
#define RESPONSEOPERATOR_HH_

#include "SP.hh"

namespace erme
{

/*!
 *  \class ResponseOperator
 *  \brief Base class for response operators
 */
class ResponseOperator
{

public:

  ResponseOperator();

  /// Virtual Destructor
  virtual ResponseOperator(){}

private:

};

} // end namespace erme

#endif // RESPONSEOPERATOR_HH_ 

//---------------------------------------------------------------------------//
//              end of file ResponseOperator.hh
//---------------------------------------------------------------------------//
