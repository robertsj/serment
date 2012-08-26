//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Node.hh
 * \brief  Node class definition
 * \author Jeremy Roberts
 * \date   Aug 22, 2012
 */
//---------------------------------------------------------------------------//

#ifndef NODE_HH_
#define NODE_HH_

#include "serment_config.h"

// ERME geometry
#include "Point.hh"

// Serment Comm
#include "Comm.hh"

// Detran
#include "DBC.hh"
#include "SP.hh"

// System
#include <vector>
#ifdef SERMENT_ENABLE_BOOST
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>
#endif

namespace erme_geometry
{

/*!
 *  \class Node
 *  \brief Base node class
 *
 *  The response matrix method decomposes a system into independent
 *  computational nodes.  This class represents an abstract node,
 *  from which nodes with specific geometries can be derived.
 *
 *  A Node represents everything needed by the underlying
 *  response generator, be it Detran or something else.  Note,
 *  the maximum response order per variable is specified.  However,
 *  the actual order realized could depend on dynamic order
 *  schemes.
 *
 *  Geometrically, nodes are rather flexible.  In 1D, they can
 *  represent slab slices, cylindrical shells, or spherical shells.
 *  In 2D, any closed shape can be well-defined.  In 3D, arbitrary
 *  right cylinders can be defined.  For 3D, south and north surfaces
 *  (in the xy plane) are enforced to simplify polarity shifts
 *  in the polar angle upon reflection at global boundaries.  The
 *  south and north surface must be the last two surfaces, respectively,
 *  in any 3D node.
 *
 */
class Node
{

public:

  /// Boundary condition identifiers
  const static int REFLECT = -1;
  const static int VACUUM  = -2;

  /// Useful typedefs
  typedef detran::SP<Node>            SP_node;
  typedef unsigned int                size_type;
  typedef std::vector<int>            vec_int;
  typedef std::vector<double>         vec_dbl;
  typedef std::vector<size_type>      vec_size_type;
  typedef std::vector<vec_size_type>  vec2_size_type;
  typedef detran::Point               Point;

  /*!
   *  \brief Constructor
   */
  Node(const size_type  dim,
       const size_type  num_surface,
       const int        nodeid,
       std::string      nodename,
       const Point      nodeorigin,
       vec2_size_type   so,
       vec_size_type    po,
       vec_size_type    ao,
       vec_size_type    eo);

  /// Virtual destructor
  virtual ~Node(){}

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Return the area of a node surface
   *  \param surface    Surface index
   */
  virtual double area(const size_type surface) const = 0;

  /// Return the volume of the node
  virtual double volume() const = 0;

  /// Return the color associated with the spatial coordinate.
  virtual double color(Point point) = 0;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  size_type dimension() const
  {
    return d_dimension;
  }

  size_type number_surfaces() const
  {
    return d_number_surfaces;
  }

  int id() const
  {
    return d_id;
  }

  std::string name() const
  {
    return d_name;
  }

  Point origin() const
  {
    return d_origin;
  }

  size_type spatial_order(const size_type s, const size_type d) const
  {
    // Preconditions
    Require(s < d_number_surfaces);
    Require(d < d_spatial_order[s].size());
    return d_spatial_order[s][d];
  }

  size_type polar_order(const size_type s) const
  {
    // Preconditions
    Require(s < d_number_surfaces);
    return d_polar_order[s];
  }

  size_type azimuthal_order(const size_type s) const
  {
    // Preconditions
    Require(s < d_number_surfaces);
    return d_azimuthal_order[s];
  }

  size_type energy_order(const size_type s) const
  {
    // Preconditions
    Require(s < d_number_surfaces);
    return d_energy_order[s];
  }

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Problem dimension
  size_type d_dimension;

  /// Number of node surfaces
  size_type d_number_surfaces;

  /// Identifier
  int d_id;

  /// Node name
  std::string d_name;

  /// Origin
  Point d_origin;

  /// Spatial orders per dimension per surface
  vec2_size_type d_spatial_order;

  /// Polar order per surface
  vec_size_type d_polar_order;

  /// Azimuth order per surface
  vec_size_type d_azimuthal_order;

  /// Energy order per surface
  vec_size_type d_energy_order;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

protected:

  /// Default constructor is needed for serialization.
  Node()
    : d_dimension(0),
      d_number_surfaces(0),
      d_id(0),
      d_origin(Point(0, 0, 0))
  {/* ... */}

private:

  friend class boost::serialization::access;

  template <class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & d_dimension;
    ar & d_number_surfaces;
    ar & d_id;
    ar & d_name;
    ar & d_spatial_order;
    ar & d_polar_order;
    ar & d_azimuthal_order;
    ar & d_energy_order;
  }

};

} // end namespace erme_geometry

BOOST_SERIALIZATION_ASSUME_ABSTRACT(erme_geometry::Node)

#endif // NODE_HH_ 

//---------------------------------------------------------------------------//
//              end of file Node.hh
//---------------------------------------------------------------------------//
