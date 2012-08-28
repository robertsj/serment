//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   NodeList.i.hh
 * \brief  NodeList inline member definitions
 * \author Jeremy Roberts
 * \date   Aug 27, 2012
 */
//---------------------------------------------------------------------------//

#ifndef NODELIST_I_HH_
#define NODELIST_I_HH_

namespace erme_geometry
{

inline void NodeList::add_node(SP_node n, vec_neighbor neighbors)
{
  Require(!is_finalized());
  Require(n);
  Require(neighbors.size() == n->number_surfaces());
  d_nodes.push_back(n);
  d_neighbors.push_back(neighbors);
}

inline NodeList::SP_node NodeList::node(const int n)
{
  Require(is_finalized());
  return d_nodes[n];
}

inline NodeList::size_type NodeList::lower_bound() const
{
  return d_lower_bound;
}

inline NodeList::size_type NodeList::upper_bound() const
{
  return d_upper_bound;
}

inline NodeList::size_type NodeList::number_global_nodes() const
{
  return d_nodes.size();
}

inline NodeList::size_type NodeList::number_local_nodes() const
{
  return d_upper_bound - d_lower_bound;
}

inline const NeighborSurface&
NodeList::neighbor(const size_type n, const size_type s) const
{
  Require(s < d_nodes[n]->number_surfaces());
  return d_neighbors[n][s];
}

// L-to-G
inline NodeList::size_type NodeList::global_index(const size_type li) const
{
  Require(li < d_upper_bound);
  size_type gi = li + d_lower_bound;
  Ensure(gi < d_nodes.size());
  return gi;
}

// G-to-L
inline int NodeList::local_index(const size_type gi) const
{
  Require(gi < d_nodes.size());
  int li = gi - d_lower_bound;
  if (li >= 0)
    return li;
  else
    return -1;
}

} // end namespace erme_geometry

#endif // NODELIST_I_HH_ 

//---------------------------------------------------------------------------//
//              end of file NodeList.i.hh
//---------------------------------------------------------------------------//
