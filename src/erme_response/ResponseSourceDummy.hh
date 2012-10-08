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
#include "ResponseSource.hh"
#include "erme_geometry/DummyNode.hh"

namespace erme_response
{

/*!
 *  \class ResponseSourceDummy
 *  \brief Fake response source for testing purposes.
 */
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

  void compute(SP_response response, ResponseIndex index)
  {
    Require(response);

    size_t in = index.nodal;

    // Easy value to recreate.
    double value = 1000000.0 * index.node +
                    100000.0 * index.surface +
                     10000.0 * index.polar +
                      1000.0 * index.azimuth +
                       100.0 * index.space0 +
                        10.0 * index.space1 +
                         1.0 * index.energy;

    for (int out = 0; out < response->size(); out++)
    {
      response->boundary_response(out, in) = value + 0.1;
    }
    response->fission_response(in) = value + 0.2;
    response->absorption_response(in) = value + 0.3;
    for (int s = 0; s < response->number_surfaces(); s++)
    {
      response->leakage_response(s, in) = value + 0.4 + 0.01 * s;
    }

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

#endif /* RESPONSESOURCEDUMMY_HH_ */
