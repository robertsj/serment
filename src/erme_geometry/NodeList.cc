//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   NodeList.cc
 *  @brief  NodeList member definitions
 *  @author Jeremy Roberts
 *  @date   Aug 27, 2012
 */
//---------------------------------------------------------------------------//

#include "NodeList.hh"
#include "NodeSerialization.hh"
#include <algorithm>

namespace erme_geometry
{

//---------------------------------------------------------------------------//
NodeList::NodeList()
  : d_lower_bound(0),
    d_upper_bound(0),
    d_number_local_surfaces(0),
    d_number_global_surfaces(0),
    d_is_finalized(false)
{
  /* ... */
}

//---------------------------------------------------------------------------//
NodeList::SP_nodelist NodeList::Create()
{
  SP_nodelist p(new NodeList());
  return p;
}

//---------------------------------------------------------------------------//
void NodeList::set_bounds(const size_t lb, const size_t ub)
{
  d_lower_bound = lb;
  d_upper_bound = ub;

  // Get the unique elements in the node map between lb and ub.  First,
  // sort, and then use unique to get the local unique nodes that actually
  // have responses.
  d_unique_nodes.resize(ub-lb);
  for (int i = lb; i < lb; ++i) d_unique_nodes[i-lb] = d_node_map[i];
  std::sort(d_unique_nodes.begin(), d_unique_nodes.end());
  vec_int::iterator it = std::unique(d_unique_nodes.begin(), d_unique_nodes.end());
  d_unique_nodes.erase(it, d_unique_nodes.end());

}

//---------------------------------------------------------------------------//
void NodeList::add_node(SP_node node)
{
  // Preconditions
  Require(!is_finalized());
  Require(node);

  d_nodes.push_back(node);
}

//---------------------------------------------------------------------------//
void NodeList::set_nodal_map(const vec_int &node_map,
                             const vec2_neighbor &neighbors)
{
  // Preconditions
  Require(is_finalized());
  Require(node_map.size() == neighbors.size());
  for (int i = 0; i < node_map.size(); ++i)
  {
    Require(node_map[i] < d_nodes.size());
    Require(d_nodes[node_map[i]]);
    Require(neighbors[i].size() == d_nodes[node_map[i]]->number_surfaces());
  }

  d_node_map = node_map;
  d_neighbors = neighbors;
}

} // end namespace erme_geometry

//---------------------------------------------------------------------------//
//              end of file NodeList.cc
//---------------------------------------------------------------------------//
