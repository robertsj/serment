//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseIndexer.i.hh
 * \brief  ResponseIndexer inline member definitions
 * \author Jeremy Roberts
 * \date   Aug 24, 2012
 */
//---------------------------------------------------------------------------//

#ifndef RESPONSEINDEXER_I_HH_
#define RESPONSEINDEXER_I_HH_

namespace erme_response
{

inline ResponseIndexer::size_t
ResponseIndexer::number_nodes() const
{
  return d_sizes.size();
}

inline ResponseIndexer::size_t
ResponseIndexer::number_surface_moments(const size_t node,
                                        const size_t surface) const
{
  Require(node < d_indices.size());
  Require(surface < d_indices[node].size());
  return d_indices[node][surface].size();
}

inline ResponseIndexer::size_t
ResponseIndexer::number_node_moments(const size_t global_node) const
{
  Require(global_node < d_sizes.size());
  return d_sizes[global_node];
}

inline ResponseIndexer::size_t
ResponseIndexer::number_local_moments() const
{
  return d_local_size;
}

inline ResponseIndexer::size_t
ResponseIndexer::number_global_moments() const
{
  return d_global_size;
}

inline ResponseIndex
ResponseIndexer::node_index(const size_t node,
                            const size_t surface,
                            const size_t index) const
{
  Require(node    < d_indices.size());
  Require(surface < d_indices[node].size());
  Require(index   < d_indices[node][surface].size());
  return d_indices[node][surface][index];
}

inline ResponseIndex
ResponseIndexer::node_index(const size_t index) const
{
  Require(index < d_local_size);
  size_t node    = d_local_indices[index][0];
  size_t surface = d_local_indices[index][1];
  size_t nindex  = d_local_indices[index][2];
  return d_indices[node][surface][nindex];
}

} // end namespace erme_response

#endif // RESPONSEINDEXER_I_HH_ 

//---------------------------------------------------------------------------//
//              end of file ResponseIndexer.i.hh
//---------------------------------------------------------------------------//
