//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   NodeList.hh
 *  @brief  NodeList class definition
 *  @author Jeremy Roberts
 *  @date   Aug 22, 2012
 */
//---------------------------------------------------------------------------//

#ifndef erme_geometry_NODELIST_HH_
#define erme_geometry_NODELIST_HH_

#include "Node.hh"
#include "NeighborSurface.hh"
#include "utilities/SP.hh"
#include "utilities/Point.hh"
#include <vector>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/access.hpp>

namespace erme_geometry
{

// Forward declare the partitioner
class NodePartitioner;

/**
 *  @class NodeList
 *  @brief Container for nodes and indices of their neighbors.
 *
 *  The node list maintains a vector of unique nodes, their
 *  placement in the global domain, and their neighbors.  For
 *  now, a "unique" node is defined by both its underlying
 *  transport problem *and* its requested expansion orders.
 *  Hence, if the same node is used multiple times in a problem
 *  with varying orders (e.g. for scoping adaptivity), then
 *  those are unique.
 *
 *  A neighbor index must be supplied for each nodal surface.
 *  If the surface is a global
 *  boundary, then either Node::REFLECT or Node::VACUUM must
 *  be specified.  The surfaces are in a Node-specific
 *  ordering.  For example, Cartesian nodes are indexed
 *  as follows: left, right, bottom, top, south, north.
 *
 *  Once the user adds all nodes, finalize() must be
 *  called.
 */
/**
 *  @example erme_geometry/test/test_NodeList.cc
 *
 *  Example of the NodeList class.
 */

class NodeList
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<NodeList>    SP_nodelist;
  typedef unsigned int                      size_t;
  typedef Node::SP_node                     SP_node;
  typedef std::vector<SP_node>              vec_node;
  typedef std::vector<int>                  vec_int;
  typedef std::vector<vec_int>              vec2_int;
  typedef std::vector<NeighborSurface>      vec_neighbor;
  typedef std::vector<vec_neighbor>         vec2_neighbor;
  typedef detran_utilities::Point           Point;
  typedef std::vector<Point>                vec_point;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /// Constructor
  NodeList();

  /// SP constructor
  static SP_nodelist Create();

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /**
   *  @brief Set the local node array bounds
   *
   *  The entire vector of nodes lives on all processes.  The
   *  partitioner must set the bounds corresponding to a
   *  particular node.  The unique nodes on each process are also
   *  identified.
   *
   *  @param lb   Lower bound
   *  @param ub   Upper bound
   */
  void set_bounds(const size_t lb, const size_t ub);

  /**
   *  @brief Add a unique node
   *  @param n          Node to be added
   */
  void add_node(SP_node n);

  /**
   *  @brief Set the map of nodes with their neighbors
   *  @param nodal_indices  Indices into the vector of unique nodes
   *  @param neighbor       Neighbor data for each node
   *  @param origins        Node origins
   */
  void set_nodal_map(const vec_int &nodes,
                     const vec2_neighbor &neighbors,
                     const vec_point &origins = vec_point(0));

  /**
   *  @brief Get a node
   *  @param n  Index into vector of unique nodes via the global index
   */
  SP_node node(const int n) const;

  /**
   *  @brief Get a node
   *  @param n  Index into vector of unique nodes via the cardinal index
   */
  SP_node unique_node(const int n) const;

  /**
   *  @brief Get a node.  Returns null pointer if not found.
   *  @param name  Name of node
   */
  //SP_node node(std::string name);

  /// Return local lower bound
  size_t lower_bound() const;
  /// Return local upper bound
  size_t upper_bound() const;

  /// Number of global nodes
  size_t number_global_nodes() const;
  /// Number of local nodes
  size_t number_local_nodes() const;
  /// Number of local surfaces
  size_t number_global_surfaces() const;
  /// Number of local surfaces
  size_t number_local_surfaces() const;
  /// Number of unique global nodes
  size_t number_unique_global_nodes() const;
  /// Number of unique local nodes
  size_t number_unique_local_nodes() const;

  /**
   *  @brief Get the global neighbor index for a node surface
   *  @param node_g  Global node index
   *  @param s  Node surface
   */
  const NeighborSurface& neighbor(const size_t node_g, const size_t s) const;

  /**
   *  @brief Get the global index of a local node
   *
   *  This indexes the nodal map that defines the geometry
   *  of the problem.
   *
   *  @param node_l  Local node index
   */
  size_t global_index_from_local(const size_t node_l) const;

  /**
   *  @brief Get the local index of a global node
   *
   *  Global index translated to the local portion of
   *  the nodal map. Returns a negative value if the
   *  global index is not within the local range.
   *
   *  @param node_g  Global node index
   */
  int local_index_from_global(const size_t node_g) const;

  /**
   *  @brief Get the unique index of a global node
   *
   *  This simply returns the nodal map entry.  Hence, the
   *  entries in the map must correspond to the nodes actually
   *  defined by the user.
   *
   *  @param n  Global node index
   */
  size_t unique_global_index_from_global(const size_t node_g) const;

  /**
   *  @brief Get the unique local index from the unique global index
   *
   *  Starting with the complete nodal map, the portion in this
   *  local group is sorted and unique elements found.  This
   *  indexer searches for the location within that sorted
   *  segment for the given global identifier.  If not found,
   *  returns negative.
   *
   *  @param node_ug  Unique global node index
   */
  int unique_local_index_from_unique_global(const size_t node_ug) const;

  /**
   *  @brief Get the unique global index from the unique local index
   *  @param node_ul  Unique local index
   */
  size_t unique_global_index_from_unique_local(const size_t node_ul) const;

  /// Have all nodes been added?
  bool is_finalized() const
  {
    return d_is_finalized;
  }

  /// Display all the nodes in the list
  void display() const;

  /// Let partitioner finalize the node list
  friend class NodePartitioner;

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Vector of unique nodes
  vec_node d_nodes;
  /// Global node map that provides indices into the node vector
  vec_int d_node_map;
  /// Indices of unique nodes needed locally
  vec_int d_unique_nodes;
  /// Vector of vectors of node (neighbor, surface) index pairs
  vec2_neighbor d_neighbors;
  /// Vector of node origins
  vec_point d_origins;
  /// Starting index of local nodes
  size_t d_lower_bound;
  /// Upper bound of local nodes
  size_t d_upper_bound;
  /// Number of local surfaces
  size_t d_number_local_surfaces;
  /// Number of global surfaces
  size_t d_number_global_surfaces;
  /// Flag indicating ready to use
  bool d_is_finalized;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /// Signify that all nodes are added.  (The partitioner calls this)
  void finalize()
  {
    d_is_finalized = true;
  }

  //-------------------------------------------------------------------------//
  // SERIALIZATION
  //-------------------------------------------------------------------------//

  friend class boost::serialization::access;

  template <class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & d_nodes;
    ar & d_node_map;
    ar & d_unique_nodes;
    ar & d_neighbors;
    ar & d_lower_bound;
    ar & d_upper_bound;
    ar & d_is_finalized;
  }

};

} // end namespace erme_geometry

//-------------------------------------------------------------------------//
// INLINE MEMBER DEFINITIONS
//-------------------------------------------------------------------------//

#include "NodeList.i.hh"

#endif // erme_geometry_NODELIST_HH_

//---------------------------------------------------------------------------//
//              end of file NodeList.hh
//---------------------------------------------------------------------------//
