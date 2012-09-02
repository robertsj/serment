//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DummyNode.cc
 * \brief  DummyNode member definitions
 * \author Jeremy Roberts
 * \date   Sep 1, 2012
 */
//---------------------------------------------------------------------------//

#include "DummyNode.hh"

namespace erme_geometry
{

CartesianNodeDummy::CartesianNodeDummy(const size_type dim,
                                       const size_type id,
                                       const size_type so,
                                       const size_type po,
                                       const size_type ao,
                                       const size_type eo)
  : CartesianNode(dim,
                  2*dim,
                  id,
                  "dummy_cartesian",
                  Point(0, 0, 0),
                  vec2_size_type(2*dim, vec_size_type(dim-1, so)),
                  vec_size_type(2*dim, po),
                  vec_size_type(2*dim, ao),
                  vec_size_type(2*dim, eo),
                  vec_dbl(3, 1.0))
{
  /* ... */
}


} // end namespace erme_geometry

#ifdef SERMENT_ENABLE_BOOST
BOOST_CLASS_EXPORT_IMPLEMENT(erme_geometry::CartesianNodeDummy)
#endif

//---------------------------------------------------------------------------//
//              end of file DummyNode.cc
//---------------------------------------------------------------------------//
