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

  typedef CartesianNode Base;

  CartesianNodeDummy(const size_type dim,
                     const size_type id = 0,
                     const size_type so = 0,
                     const size_type po = 0,
                     const size_type ao = 0,
                     const size_type eo = 0);

  double color(Point point)
  {
    return 0.0;
  }

private:

  /// Default constructor needed for serialization
  CartesianNodeDummy(){}

#ifdef SERMENT_ENABLE_BOOST
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Base>(*this);
  }
#endif

};

} // end namespace erme_geometry

#ifdef SERMENT_ENABLE_BOOST
BOOST_CLASS_EXPORT_KEY(erme_geometry::CartesianNodeDummy)
#endif

#endif // DUMMYNODE_HH_ 

//---------------------------------------------------------------------------//
//              end of file DummyNode.hh
//---------------------------------------------------------------------------//
