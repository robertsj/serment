//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  HexagonalNode.hh
 *  @brief HexagonalNode class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef erme_geometry_HEXAGONALNODE_HH_
#define erme_geometry_HEXAGONALNODE_HH_

#include "Node.hh"

namespace erme_geometry
{

/**
 *  @class HexagonalNode
 *  @brief Base Hexagonal node class
 *
 */
class HexagonalNode: public Node
{

public:

  //--------------------------------------------------------------------------//
  // ENUMERATIONS
  //--------------------------------------------------------------------------//

  // Face identifiers
  enum FACE
  {
    SOUTHWEST,
    NORTHEAST,
    NORTHWEST,
    SOUTHEAST,
    SOUTH,
    NORTH,
    BOTTOM,
    TOP
  };

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef Node                    Base;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  HexagonalNode(const size_t  dimension,
                std::string   nodename = "hexagonal_node",
                vec2_size_t   so,
                vec_size_t    po,
                vec_size_t    ao,
                vec_size_t    eo,
                const double  face_width,
                const double  height = 1.0);

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE
  //--------------------------------------------------------------------------//

  /**
   *  @brief Return the area of a node surface
   *  @param surface    Surface index
   */
  double area(const size_t surface) const
  {
    Require(surface < number_surfaces());
    if (surface < 6)
      return d_face_width * d_height;
    // area of xy faces = 3*sqrt(3)/2 * face_width^2
    return 2.59807621135332 * d_face_width * d_face_width;
  }

  /// Return the volume of the node
  double volume() const
  {
    return 2.59807621135332 * d_face_width * d_face_width * d_height;
  }

  /// Return the color associated with the spatial coordinate.
  virtual double color(const Point &point, std::string key = "MATERIAL") = 0;

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Width of a face
  const double d_face_width;
  /// Height of cell (defaults to 1.0 for 2D)
  const double d_height;


  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

protected:

  /// Default constructor needed for serialization
  HexagonalNode() : d_face_width(0.0), d_height(0.0) {}

private:

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Base>(*this);
    ar & d_face_width;
    ar & d_height;
  }

};

} // end namespace erme_geometry

//BOOST_CLASS_EXPORT_KEY(erme_geometry::HexagonalNode)
//BOOST_SERIALIZATION_ASSUME_ABSTRACT(erme_geometry::HexagonalNode)

#endif // erme_geometry_HEXAGONALNODE_HH_

//----------------------------------------------------------------------------//
//              end of file HexagonalNode.hh
//----------------------------------------------------------------------------//
