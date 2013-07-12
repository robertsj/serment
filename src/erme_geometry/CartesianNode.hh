//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  CartesianNode.hh
 *  @brief CartesianNode class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef erme_geometry_CARTESIANNODE_HH_
#define erme_geometry_CARTESIANNODE_HH_

#include "Node.hh"

namespace erme_geometry
{

/**
 *  @class CartesianNode
 *  @brief Base Cartesian node class
 */
class CartesianNode: public Node
{

public:

  //--------------------------------------------------------------------------//
  // ENUMERATIONS
  //--------------------------------------------------------------------------//

  // Face identifiers
  enum FACE
  {
    WEST, EAST, SOUTH, NORTH, BOTTOM, TOP, END_FACE
  };

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef Node                    Base;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  CartesianNode(const size_t  dimension,
                std::string   nodename,
                vec2_size_t   so,
                vec_size_t    po,
                vec_size_t    ao,
                vec_size_t    eo,
                vec_dbl       nodewidth);

  /// SP constructor
  static SP_node Create(const size_t  dimension,
                        std::string   nodename,
                        vec2_size_t   so,
                        vec_size_t    po,
                        vec_size_t    ao,
                        vec_size_t    eo,
                        vec_dbl       nodewidth);

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE
  //--------------------------------------------------------------------------//

  /**
   *  @brief Return the area of a node surface
   *  @param surface    Surface index
   */
  double area(const size_t surface) const;

  /// Return the volume of the node
  double volume() const;

  /// Default color.  Am I in the box or not?
  double color(const Point &point, std::string = "MATERIAL");

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  double width(const size_t dim) const;

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Width in each dimension
  vec_dbl d_width;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

protected:

  /// Default constructor needed for serialization
  CartesianNode(){}

private:

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Base>(*this);
    ar & d_width;
  }

};

} // end namespace erme_geometry

#endif // erme_geometry_CARTESIANNODE_HH_

//----------------------------------------------------------------------------//
//              end of file CartesianNode.hh
//----------------------------------------------------------------------------//
