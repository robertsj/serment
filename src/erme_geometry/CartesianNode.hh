//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   CartesianNode.hh
 *  @brief  CartesianNode
 *  @author Jeremy Roberts
 *  @date   Aug 22, 2012
 */
//---------------------------------------------------------------------------//

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

  //-------------------------------------------------------------------------//
  // ENUMERATIONS
  //-------------------------------------------------------------------------//

  // Face identifiers
  enum FACE
  {
    WEST,
    EAST,
    SOUTH,
    NORTH,
    BOTTOM,
    TOP
  };

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef Node                    Base;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  CartesianNode(const size_t  d,
                const size_t  n,
                const int     nodeid,
                std::string   nodename,
                const Point   nodeorigin,
                vec2_size_t   so,
                vec_size_t    po,
                vec_size_t    ao,
                vec_size_t    eo,
                vec_dbl       nodewidth);

  /// SP constructor
  static SP_node Create(const size_t  d,
                        const size_t  n,
                        const int     nodeid,
                        std::string   nodename,
                        const Point   nodeorigin,
                        vec2_size_t   so,
                        vec_size_t    po,
                        vec_size_t    ao,
                        vec_size_t    eo,
                        vec_dbl       nodewidth)
  {
    SP_node p(new CartesianNode(d, n, nodeid, nodename, nodeorigin,
                                so, po, ao, eo, nodewidth));
    return p;
  }


  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE
  //-------------------------------------------------------------------------//

  /**
   *  @brief Return the area of a node surface
   *  @param surface    Surface index
   */
  double area(const size_t surface) const
  {
    Require(surface < number_surfaces());
    if (surface == 0 or surface == 1)
      return d_width[1] * d_width[2];
    if (surface == 2 or surface == 3)
      return d_width[0] * d_width[2];
    return d_width[0] * d_width[1];
  }

  /// Return the volume of the node
  double volume() const
  {
    return d_width[0] * d_width[1] * d_width[2];
  }

  /// Default color.  Am I in the box or not?
  virtual double color(Point point)
  {
    Point point_local = point + -1*origin();
    if ( (point_local.x() < 0.0) or
         (point_local.y() < 0.0) or
         (point_local.z() < 0.0) or
         (point_local.x() > d_width[0]) or
         (point_local.y() > d_width[1]) or
         (point_local.z() > d_width[2]) )
    {
      return -1.0;
    }
    return (double) id();
  }


private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Width in each dimension
  vec_dbl d_width;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

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

//---------------------------------------------------------------------------//
//              end of file CartesianNode.hh
//---------------------------------------------------------------------------//
