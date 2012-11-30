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

/**
 *  Dummy Cartesian node for testing
 */
class CartesianNodeDummy: public CartesianNode
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef CartesianNode Base;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  CartesianNodeDummy(const size_t dim,
                     const size_t so = 0,
                     const size_t po = 0,
                     const size_t ao = 0,
                     const size_t eo = 0);

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE
  //-------------------------------------------------------------------------//

  double color(Point point)
  {
    return 0.0;
  }

private:

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /// Default constructor needed for serialization
  CartesianNodeDummy(){}

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Base>(*this);
  }

};

} // end namespace erme_geometry

#endif // DUMMYNODE_HH_ 

//---------------------------------------------------------------------------//
//              end of file DummyNode.hh
//---------------------------------------------------------------------------//
