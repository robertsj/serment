//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   NodeList.hh
 * \brief  NodeList class definition
 * \author Jeremy Roberts
 * \date   Aug 22, 2012
 */
//---------------------------------------------------------------------------//

#ifndef NODELIST_HH_
#define NODELIST_HH_

#include "Node.hh"
#include <vector>

namespace erme_geometry
{

/*!
 *  \class NodeList
 *  \brief Container for nodes and indices of their neighbors.
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

class NodeList
{

public:

  typedef unsigned int          size_type;
  typedef Node::SP_node         SP_node;
  typedef std::vector<SP_node>  vec_node;
  typedef std::vector<int>      vec_int;
  typedef std::vector<vec_int>  vec2_int;

  /// Constructor
  NodeList() : d_is_finalized(false), d_global_offset(0) {}

  void add_node(SP_node n, vec_int neighbors)
  {
    Require(!is_finalized());
    Require(n);
    Require(neighbors.size() == n->number_surfaces());
    d_nodes.push_back(n);
    d_neighbors.push_back(neighbors);
  }

  void resize(const size_type n)
  {
    d_nodes.resize(n);
    d_neighbors.resize(n);
    d_is_finalized = false;
  }

  SP_node node(const int n)
  {
    Require(is_finalized());
    Require(n < number_nodes());
    return d_nodes[n];
  }

  size_type number_nodes() const
  {
    Require(is_finalized());
    return d_nodes.size();
  }

  /*!
   *  \brief Get the neighbor index for a node surface
   *  \param n  Node index
   *  \param s  Node surface
   */
  int neighbor(const size_type n, const size_type s)
  {
    Require(is_finalized());
    Require(n < number_nodes());
    Require(s < d_nodes[n]->number_surfaces());
    return d_neighbors[n][s];
  }

  vec_int& neighbor(const size_type n)
  {
    Require(is_finalized());
    Require(n < number_nodes());
    return d_neighbors[n];
  }

  /*!
   *  \brief Get the global index of a node
   *  \param n local node index
   */
  size_type global_index(const size_type n)
  {
    return n + d_global_offset;
  }

  void finalize()
  {
    d_is_finalized = true;
  }

  bool is_finalized() const
  {
    return d_is_finalized;
  }

private:

  /// List of nodes
  vec_node d_nodes;

  /// List of node neighbor vectors
  vec2_int d_neighbors;

  /// Finalized?
  bool d_is_finalized;

  /// Difference between local node index and its global index
  size_type d_global_offset;

};

} // end namespace erme_geometry

#endif // NODELIST_HH_ 

//---------------------------------------------------------------------------//
//              end of file NodeList.hh
//---------------------------------------------------------------------------//
