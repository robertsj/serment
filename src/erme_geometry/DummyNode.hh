//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DummyNode.hh
 * \brief  DummyNode 
 * \author Jeremy Roberts
 * \date   Aug 29, 2012
 */
//---------------------------------------------------------------------------//

#ifndef DUMMYNODE_HH_
#define DUMMYNODE_HH_

#include "erme_geometry/CartesianNode.hh"

namespace erme_geometry
{

/*!
 *  Dummy Cartesian node for testing
 */
class CartesianNodeDummy: public CartesianNode
{

public:

  CartesianNodeDummy(const size_type dim,
                     const size_type so = 0,
                     const size_type po = 0,
                     const size_type ao = 0,
                     const size_type eo = 0)
    : CartesianNode(dim,
                    2*dim,
                    1234,
                    "dummy_cartesian",
                    Point(0, 0, 0),
                    vec2_size_type(2*dim, vec_size_type(dim-1, so)),
                    vec_size_type(2*dim, po),
                    vec_size_type(2*dim, ao),
                    vec_size_type(2*dim, eo),
                    vec_dbl(dim, 1.0))
  {
    /* ... */
  }

  double color(Point ppint)
  {
    return 0.0;
  }

};

} // end namespace erme_geometry

#endif // DUMMYNODE_HH_ 

//---------------------------------------------------------------------------//
//              end of file DummyNode.hh
//---------------------------------------------------------------------------//
