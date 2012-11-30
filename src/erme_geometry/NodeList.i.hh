//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   NodeList.i.hh
 *  @brief  NodeList inline member definitions
 *  @author Jeremy Roberts
 *  @date   Aug 27, 2012
 */
//---------------------------------------------------------------------------//

#ifndef erme_geometry_NODELIST_I_HH_
#define erme_geometry_NODELIST_I_HH_

namespace erme_geometry
{

//---------------------------------------------------------------------------//
inline NodeList::SP_node NodeList::node(const int n) const
{
  // Preconditions
  Require(is_finalized());
  Require(n < number_global_nodes());

  return d_nodes[d_node_map[n]];
}

//---------------------------------------------------------------------------//
inline NodeList::size_t NodeList::lower_bound() const
{
  return d_lower_bound;
}

//---------------------------------------------------------------------------//
inline NodeList::size_t NodeList::upper_bound() const
{
  return d_upper_bound;
}

//---------------------------------------------------------------------------//
inline NodeList::size_t NodeList::number_global_nodes() const
{
  return d_node_map.size();
}

//---------------------------------------------------------------------------//
inline NodeList::size_t NodeList::number_local_nodes() const
{
  return d_upper_bound - d_lower_bound;
}

//---------------------------------------------------------------------------//
inline NodeList::size_t NodeList::number_global_surfaces() const
{
  size_t ns = 0;
  for (int n = 0; n < number_global_nodes(); n++)
    ns += node(n)->number_surfaces();
  return ns;
}

//---------------------------------------------------------------------------//
inline NodeList::size_t NodeList::number_local_surfaces() const
{
  size_t ns = 0;
  for (int n = lower_bound(); n < upper_bound(); n++)
    ns += node(n)->number_surfaces();
  return ns;
}

//---------------------------------------------------------------------------//
inline const NeighborSurface&
NodeList::neighbor(const size_t n, const size_t s) const
{
  // Preconditions
  Require(s < node(n)->number_surfaces());

  return d_neighbors[n][s];
}

//---------------------------------------------------------------------------//
// L-to-G
inline NodeList::size_t NodeList::global_index(const size_t li) const
{
  // Preconditions
  Require(li < d_upper_bound);

  size_t gi = li + d_lower_bound;

  // Postconditions
  Ensure(gi < number_global_nodes());
  return gi;
}

//---------------------------------------------------------------------------//
// G-to-L
inline int NodeList::local_index(const size_t gi) const
{
  // Preconditions
  Require(gi < d_nodes.size());

  int li = gi - d_lower_bound;
  if (li >= 0)
    return li;
  else
    return -1;
}

} // end namespace erme_geometry

#endif // erme_geometry_NODELIST_I_HH_

//---------------------------------------------------------------------------//
//              end of file NodeList.i.hh
//---------------------------------------------------------------------------//
