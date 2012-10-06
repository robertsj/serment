//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   NodeList.hh
 *  @brief  NodeList class definition
 *  @author Jeremy Roberts
 *  @date   Aug 22, 2012
 */
//---------------------------------------------------------------------------//

#ifndef NODELIST_HH_
#define NODELIST_HH_

#include "Node.hh"
#include "NeighborSurface.hh"
#include "utilities/SP.hh"
#include <vector>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/access.hpp>

namespace erme_geometry
{

/**
 *  @class NodeList
 *  @brief Container for nodes and indices of their neighbors.
 *
 *  The purpose of NodeList is to contain a list of Node
 *  pointers and their neighbors.  A neighbor index must be
 *  supplied for each surface.  If the surface is a global
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

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /// Constructor
  NodeList();

  /// SP constructor
  static SP_nodelist Create();

  /**
   *  @brief Set the local node array bounds
   *
   *  The entire vector of nodes lives on all process.  The
   *  partitioner must set the bounds corresponding to a
   *  particular node.
   *
   *  @param lb   Lower bound
   *  @param ub   Upper bound
   */
  void set_bounds(const size_t lb, const size_t ub)
  {
    d_lower_bound = lb;
    d_upper_bound = ub;
  }

  /**
   *  @brief Add a node and a vector of (neighbor, surface) indices
   *  @param n          Node to be added
   *  @param neighbors  Vector of neighbor/surface indices
   */
  void add_node(SP_node n, vec_neighbor neighbors);

  /**
   *  @brief Get a node
   *  @param n  Local node index
   */
  SP_node node(const int n);

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

  /**
   *  @brief Get the global neighbor index for a node surface
   *  @param n  Global node index
   *  @param s  Node surface
   */
  const NeighborSurface& neighbor(const size_t n, const size_t s) const;

  /**
   *  @brief Get the global index of a local node
   *  @param n  Local node index
   */
  size_t global_index(const size_t n) const;

  /**
   *  @brief Get the local index of a global node
   *  @param n  Global node index
   *
   *  Returns a negative value if the global index
   *  is not within the local range.
   */
  int local_index(const size_t n) const;

  void finalize()
  {
    d_is_finalized = true;
  }

  bool is_finalized() const
  {
    return d_is_finalized;
  }

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Vector of nodes
  vec_node d_nodes;
  /// Vector of vectors of node (neighbor, surface) index pairs
  vec2_neighbor d_neighbors;
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
  // SERIALIZATION
  //-------------------------------------------------------------------------//

  friend class boost::serialization::access;

  template <class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & d_nodes;
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

#endif // NODELIST_HH_ 

//---------------------------------------------------------------------------//
//              end of file NodeList.hh
//---------------------------------------------------------------------------//
