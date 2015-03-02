//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Node.hh
 *  @brief Node class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef erme_geometry_NODE_HH_
#define erme_geometry_NODE_HH_

#include "serment_config.h"
#include "comm/Comm.hh"
#include "geometry/Point.hh"
#include "utilities/Definitions.hh"
#include "utilities/DBC.hh"
#include "utilities/SP.hh"
#include <vector>
#include <string>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>

/**
 *  @namespace erme_geometry
 *  @brief Contains constructs for describing problem geometry
 */
namespace erme_geometry
{

/**
 *  @class Node
 *  @brief Base node class
 *
 *  The response matrix method decomposes a system into independent
 *  computational nodes.  This class represents an abstract node,
 *  from which nodes with specific geometries can be derived.
 *
 *  A Node represents everything needed by the underlying
 *  response generator, be it Detran or something else.  Note,
 *  the maximum response order per variable is specified.  However,
 *  the actual order realized could depend on dynamic order
 *  schemes.  Hence, while a unique node can be specified once,
 *  the actual expansion orders used for multiple instances of the
 *  node may vary.
 *
 *  Geometrically, nodes are rather flexible.  In 1D, they can
 *  represent slab slices, cylindrical shells, or spherical shells.
 *  In 2D, any closed shape can be well-defined.  In 3D, arbitrary
 *  right cylinders can be defined.  For 3D, bottom and top surfaces
 *  (in the xy plane) are enforced to simplify polarity shifts
 *  in the polar angle upon reflection at global boundaries.  The
 *  bottom and top surface must be the last two surfaces, respectively,
 *  in any 3D node.
 *
 *  Because visualization of a node (or the set of connected nodes)
 *  is extremely useful, a Node can implement a \ref color function
 *  that returns some useful information about the node.  Because
 *  all nodes are assigned an origin, the implementation can determine
 *  whether a queried point is within the node, and if so, return a
 *  color.  By default, the south-west-bottom corner of the global
 *  domain is the global origin, and the s-w-b corner of nodes
 *  should be the local origin---different nodes might need to be
 *  flexible with that notion.
 *
 */
class Node
{

public:

  //--------------------------------------------------------------------------//
  // ENUMERATIONS
  //--------------------------------------------------------------------------//

  /// Boundary condition identifiers.  Note, periodic is not listed, since
  /// nodes can be explicitly linked in a periodic fashion.
  enum NODE_BC
  {
    REFLECT = -1, VACUUM = -2
  };

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<Node>        SP_node;
  typedef detran_utilities::size_t          size_t;
  typedef detran_utilities::vec_int         vec_int;
  typedef detran_utilities::vec_dbl         vec_dbl;
  typedef detran_utilities::vec_size_t      vec_size_t;
  typedef detran_utilities::vec2_size_t     vec2_size_t;
  typedef detran_geometry::Point            Point;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param dimension        Dimension of the node
   *  @param number_surfaces  Number of surfaces
   *  @param name             Unique Name
   *  @param so               Spatial orders [surface][axis]
   *  @param po               Polar orders [surface]
   *  @param ao               Azimuth orders [surface]
   *  @param eo               Energy orders [surface]
   */
  Node(const size_t  dimension,
       const size_t  number_surfaces,
       std::string   name,
       vec2_size_t   so,
       vec_size_t    po,
       vec_size_t    ao,
       vec_size_t    eo);

  /// Pure virtual destructor
  virtual ~Node() = 0;

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL NODES MUST IMPLEMENT THESE
  //--------------------------------------------------------------------------//

  /**
   *  @brief Return the area of a node surface
   *  @param surface    Surface index
   */
  virtual double area(const size_t surface) const = 0;

  /// Return the volume of the node
  virtual double volume() const = 0;

  /// Return a color
  virtual double color(const Point &point, std::string key = "MATERIAL") = 0;

  //--------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //--------------------------------------------------------------------------//

  /// Node dimension
  size_t dimension() const;
  /// Number nodal surfaces
  size_t number_surfaces() const;
  /// Number nodal pins
  size_t number_pins() const;
  /// Set number of pins
  void set_number_pins(const size_t n);
  /// Nodal name
  std::string name() const;
  /// Spatial order for a surface and possible dimension
  size_t spatial_order(const size_t s, const size_t d) const;
  /// Polar order for a surface
  size_t polar_order(const size_t s) const;
  /// Azimuthal order for a surface
  size_t azimuthal_order(const size_t s) const;
  /// Energy order for a surface
  size_t energy_order(const size_t s) const;
  /// Set fuel flag
  void set_fuel(bool flag);
  /// Is is fuel?
  bool is_fuel() const {return d_is_fuel;}
  /// Pretty print of key characteristics
  void display() const;

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Problem dimension
  size_t d_dimension;
  /// Number of node surfaces
  size_t d_number_surfaces;
  /// Identifier
  std::string d_name;
  /// Spatial orders per dimension per surface
  vec2_size_t d_spatial_order;
  /// Polar order per surface
  vec_size_t d_polar_order;
  /// Azimuth order per surface
  vec_size_t d_azimuthal_order;
  /// Energy order per surface
  vec_size_t d_energy_order;
  /// Flag indicated node is fueled or not (default: true)
  bool d_is_fuel;

protected:

  /// Number of node pins
  size_t d_number_pins;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  /// Default constructor is needed for serialization.
  Node()
    : d_dimension(0),
      d_number_surfaces(0),
      d_number_pins(0),
      d_is_fuel(true),
      d_name("no_name_given")
  {/* ... */}

private:

  friend class boost::serialization::access;

  template <class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & d_dimension;
    ar & d_number_surfaces;
    ar & d_number_pins;
    ar & d_name;
    ar & d_spatial_order;
    ar & d_polar_order;
    ar & d_azimuthal_order;
    ar & d_energy_order;
    ar & d_is_fuel;
  }

};

} // end namespace erme_geometry

#include "Node.i.hh"

#endif // erme_geometry_NODE_HH_

//----------------------------------------------------------------------------//
//              end of file Node.hh
//----------------------------------------------------------------------------//
