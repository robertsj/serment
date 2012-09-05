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

inline NodeList::size_t NodeList::lower_bound() const
{
  return d_lower_bound;
}

inline NodeList::size_t NodeList::upper_bound() const
{
  return d_upper_bound;
}

inline NodeList::size_t NodeList::number_global_nodes() const
{
  return d_nodes.size();
}

inline NodeList::size_t NodeList::number_local_nodes() const
{
  return d_upper_bound - d_lower_bound;
}

inline NodeList::size_t NodeList::number_global_surfaces() const
{
  size_t ns = 0;
  for (int n = 0; n < number_global_nodes(); n++)
    ns += d_nodes[n]->number_surfaces();
  return ns;
}

inline NodeList::size_t NodeList::number_local_surfaces() const
{
  size_t ns = 0;
  for (int n = lower_bound(); n < upper_bound(); n++)
    ns += d_nodes[n]->number_surfaces();
  return ns;
}

inline const NeighborSurface&
NodeList::neighbor(const size_t n, const size_t s) const
{
  Require(s < d_nodes[n]->number_surfaces());
  return d_neighbors[n][s];
}

// L-to-G
inline NodeList::size_t NodeList::global_index(const size_t li) const
{
  Require(li < d_upper_bound);
  size_t gi = li + d_lower_bound;
  Ensure(gi < d_nodes.size());
  return gi;
}

// G-to-L
inline int NodeList::local_index(const size_t gi) const
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
