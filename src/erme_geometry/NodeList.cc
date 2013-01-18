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
  // Preconditions
  Require(ub >= lb); // Can have more processes than nodes

  d_lower_bound = lb;
  d_upper_bound = ub;

  // Get the unique elements in the node map between lb and ub.  First,
  // sort, and then use unique to get the local unique nodes that actually
  // have responses.
  d_unique_nodes.resize(ub-lb);
  for (int i = lb; i < ub; ++i)
  {
    //std::cout << " unique[" << i - lb << "]=" << d_node_map[i] << std::endl;
    d_unique_nodes[i-lb] = d_node_map[i];
  }
  std::sort(d_unique_nodes.begin(), d_unique_nodes.end());
  vec_int::iterator it =
    std::unique(d_unique_nodes.begin(), d_unique_nodes.end());
  d_unique_nodes.erase(it, d_unique_nodes.end());

  // d_unique_nodes now has global indices of nodes this process needs.
}

//---------------------------------------------------------------------------//
void NodeList::add_node(SP_node node)
{
  // Preconditions
  Insist(!is_finalized(), "Cannot be finalized when adding a node.")
  Insist(node, "Invalid node.");

  d_nodes.push_back(node);
}

//---------------------------------------------------------------------------//
void NodeList::set_nodal_map(const vec_int          &node_map,
                             const vec2_neighbor    &neighbors,
                             const vec_point        &origins)
{
  // Preconditions
  Insist(node_map.size() == neighbors.size(),
         "Node map size is inconsistent with size of neighbors.");
  for (int i = 0; i < node_map.size(); ++i)
  {
    Insist(node_map[i] < d_nodes.size(),
           "Node index in map is larger than number of nodes added.");
    Insist(d_nodes[node_map[i]],
           "Null node indexed by map.  Ensure all nodes added are built.");
    Insist(neighbors[i].size() == d_nodes[node_map[i]]->number_surfaces(),
           "Size of node neighbor vector inconsistent with surface count.");
  }

  d_node_map = node_map;
  d_neighbors = neighbors;
  d_origins = origins;
}

//---------------------------------------------------------------------------//
NodeList::Point NodeList::origin(const size_t node_g) const
{
  Require(node_g < number_global_nodes());
  return d_origins[node_g];
}

//---------------------------------------------------------------------------//
void NodeList::display() const
{
  if (serment_comm::Comm::rank() == 0)
  {
    for (size_t n = 0; n < number_global_nodes(); ++n)
    {
      std::cout << "---------------------" << std::endl;
      std::cout << " GLOBAL NODE " << n << std::endl;
      std::cout << "---------------------" << std::endl;
      node(n)->display();
    }
  }
}

} // end namespace erme_geometry

//---------------------------------------------------------------------------//
//              end of file NodeList.cc
//---------------------------------------------------------------------------//
