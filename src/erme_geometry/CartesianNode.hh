//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CartesianNode.hh
 * \brief  CartesianNode 
 * \author Jeremy Roberts
 * \date   Aug 22, 2012
 */
//---------------------------------------------------------------------------//

#ifndef CARTESIANNODE_HH_
#define CARTESIANNODE_HH_

// ERME Geometry
#include "Node.hh"

namespace erme_geometry
{

/*!
 *  \class CartesianNode
 *  \brief Base Cartesian node class
 *
 */
class CartesianNode: public Node
{

public:

  // Face identifiers
  enum FACE
  {
    LEFT,
    RIGHT,
    BOTTOM,
    TOP,
    SOUTH,
    NORTH
  };

  typedef Node                    Base;
  typedef std::vector<double>     vec_dbl;

  CartesianNode(const size_type  d,
                const size_type  n,
                const int        nodeid,
                std::string      nodename,
                const Point      nodeorigin,
                vec2_size_type   so,
                vec_size_type    po,
                vec_size_type    ao,
                vec_size_type    eo,
                vec_dbl          nodewidth);

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Return the area of a node surface
   *  \param surface    Surface index
   */
  double area(const size_type surface) const
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

  /// Return the color associated with the spatial coordinate.
  virtual double color(Point point) = 0;

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

#ifdef SERMENT_ENABLE_BOOST
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Base>(*this);
    ar & d_width;
  }
#endif

};

} // end namespace erme_geometry

#ifdef SERMENT_ENABLE_BOOST
BOOST_CLASS_EXPORT_KEY(erme_geometry::CartesianNode)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(erme_geometry::CartesianNode)
#endif



#endif // CARTESIANNODE_HH_ 

//---------------------------------------------------------------------------//
//              end of file CartesianNode.hh
//---------------------------------------------------------------------------//
