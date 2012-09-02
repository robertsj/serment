//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseSourceDetran.hh
 * \brief  ResponseSourceDetran class definition
 * \author Jeremy Roberts
 * \date   Sep 1, 2012
 */
//---------------------------------------------------------------------------//

#ifndef RESPONSESOURCEDETRAN_HH_
#define RESPONSESOURCEDETRAN_HH_

#include "erme_geometry/CartesianNodeDetran.hh"

namespace erme_response
{

/*!
 *  \class ResponseSourceDetran
 *  \brief Compute responses using Detran
 */
class ResponseSourceDetran: public ResponseSource
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef erme_geometry::CartesianNodeDetran::SP_node SP_node;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   */
  ResponseSourceDetran(SP_node node)
    : ResponseSource(node)
  {
    THROW("implement detran source");
  }

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE
  //-------------------------------------------------------------------------//

  void compute(SP_response response, ResponseIndex index)
  {
    /* ... */
  }

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

};

} // end namespace erme_response

#endif // RESPONSESOURCEDETRAN_HH_ 

//---------------------------------------------------------------------------//
//              end of file ResponseSourceDetran.hh
//---------------------------------------------------------------------------//
