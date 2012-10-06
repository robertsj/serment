//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   DummyNode.cc
 *  @brief  DummyNode member definitions
 *  @author Jeremy Roberts
 *  @date   Sep 1, 2012
 */
//---------------------------------------------------------------------------//

#include "DummyNode.hh"
#include "NodeSerialization.hh"

namespace erme_geometry
{

CartesianNodeDummy::CartesianNodeDummy(const size_t dim,
                                       const size_t id,
                                       const size_t so,
                                       const size_t po,
                                       const size_t ao,
                                       const size_t eo)
  : CartesianNode(dim,
                  2*dim,
                  id,
                  "dummy_cartesian",
                  Point(0, 0, 0),
                  vec2_size_t(2*dim, vec_size_t(dim-1, so)),
                  vec_size_t(2*dim, po),
                  vec_size_t(2*dim, ao),
                  vec_size_t(2*dim, eo),
                  vec_dbl(3, 1.0))
{
  /* ... */
}


} // end namespace erme_geometry

BOOST_CLASS_EXPORT_IMPLEMENT(erme_geometry::CartesianNodeDummy)

//---------------------------------------------------------------------------//
//              end of file DummyNode.cc
//---------------------------------------------------------------------------//
