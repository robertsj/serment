//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseSourceDummy.hh
 * \author robertsj
 * \date   Aug 31, 2012
 * \brief  ResponseSourceDummy class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef RESPONSESOURCEDUMMY_HH_
#define RESPONSESOURCEDUMMY_HH_

#include "NodeResponse.hh"
#include "erme_geometry/DummyNode.hh"

namespace erme_response
{

class ResponseSourceDummy: public ResponseSource
{


public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef erme_geometry::CartesianNodeDummy::SP_node SP_node;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   */
  ResponseSourceDummy(SP_node node);

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE
  //-------------------------------------------------------------------------//

  /// Get a response function, i.e. the response due to one incident condition
  virtual void compute(SP_response response, ResponseIndex index);

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  size_t d_dimension;
  size

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//



};

} // end namespace erme_response

#endif /* RESPONSESOURCEDUMMY_HH_ */
