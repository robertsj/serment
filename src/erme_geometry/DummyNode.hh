//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  DummyNode.hh
 *  @brief DummyNode class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef erme_geometry_DUMMYNODE_HH_
#define erme_geometry_DUMMYNODE_HH_

#include "erme_geometry/CartesianNode.hh"

namespace erme_geometry
{

/**
 *  Dummy Cartesian node for testing
 */
class CartesianNodeDummy: public CartesianNode
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef CartesianNode Base;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  CartesianNodeDummy(const size_t dim,
                     const size_t so = 0,
                     const size_t po = 0,
                     const size_t ao = 0,
                     const size_t eo = 0);

private:

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

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

#endif // erme_geometry_DUMMYNODE_HH_

//----------------------------------------------------------------------------//
//              end of file DummyNode.hh
//----------------------------------------------------------------------------//
