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
inline NodeList::SP_node NodeList::node(const int node_g) const
{
  Require(is_finalized());
  Require(node_g < number_global_nodes());
  return d_nodes[d_node_map[node_g]];
}

//---------------------------------------------------------------------------//
inline NodeList::SP_node NodeList::unique_node(const int node_ug) const
{
  Require(is_finalized());
  Require(node_ug < number_unique_global_nodes());
  return d_nodes[node_ug];
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
inline NodeList::size_t NodeList::number_unique_global_nodes() const
{
  return d_nodes.size();
}
//---------------------------------------------------------------------------//
inline NodeList::size_t NodeList::number_local_nodes() const
{
  return d_upper_bound - d_lower_bound;
}
//---------------------------------------------------------------------------//
inline NodeList::size_t NodeList::number_unique_local_nodes() const
{
  return d_unique_nodes.size();
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
NodeList::neighbor(const size_t node_g, const size_t s) const
{
  Require(node_g < number_global_nodes());
  Require(s < node(node_g)->number_surfaces());
  return d_neighbors[node_g][s];
}

//---------------------------------------------------------------------------//
// L-to-G
inline NodeList::size_t
NodeList::global_index_from_local(const size_t node_l) const
{
  // Preconditions
  Require(node_l < d_upper_bound);

  size_t node_g = node_l + d_lower_bound;

  // Postconditions
  Ensure(node_g < number_global_nodes());
  return node_g;
}

//---------------------------------------------------------------------------//
// G-to-L
inline int NodeList::local_index_from_global(const size_t node_g) const
{
  Require(node_g < number_global_nodes());
  int node_l = node_g - d_lower_bound;
  if (node_l >= 0)
    return node_l;
  else
    return -1;
}

//---------------------------------------------------------------------------//
// G-to-GU
inline NodeList::size_t
NodeList::unique_global_index_from_global(const size_t node_g) const
{
  Require(node_g < number_global_nodes());

  size_t node_ug = d_node_map[node_g];

  Ensure(node_ug < number_unique_global_nodes());
  return node_ug;
}

//---------------------------------------------------------------------------//
// GU-to-LU
inline int
NodeList::unique_local_index_from_unique_global(const size_t node_ug) const
{
  Require(node_ug < number_unique_global_nodes());

  int node_lg = -1;
  for (size_t i = 0; i < d_unique_nodes.size(); ++i)
  {
    //std::cout << "UNIQUE = " << d_unique_nodes[i] << " ug =" << node_ug << std::endl;
    if (d_unique_nodes[i] == node_ug) node_lg = i;
  }
  return node_lg;
}

//---------------------------------------------------------------------------//
inline NodeList::size_t
NodeList::unique_global_index_from_unique_local(const size_t node_ul) const
{
  Require(node_ul < d_unique_nodes.size());

  size_t node_ug  = d_unique_nodes[node_ul];

  Ensure(node_ug < number_unique_global_nodes());
  return node_ug;
}

} // end namespace erme_geometry

#endif // erme_geometry_NODELIST_I_HH_

//---------------------------------------------------------------------------//
//              end of file NodeList.i.hh
//---------------------------------------------------------------------------//
