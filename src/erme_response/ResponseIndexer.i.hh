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

//---------------------------------------------------------------------------//

inline ResponseIndexer::size_t
ResponseIndexer::number_nodes() const
{
  return d_sizes.size();
}

inline ResponseIndexer::size_t
ResponseIndexer::number_surface_moments(const size_t node_g,
                                        const size_t surface) const
{
  Require(node_g < d_indices.size());
  Require(surface < d_indices[node_g].size());
  return d_indices[node_g][surface].size();
}

inline ResponseIndexer::size_t
ResponseIndexer::number_node_moments(const size_t node_g) const
{
  Require(node_g < d_sizes.size());
  return d_sizes[node_g];
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

//---------------------------------------------------------------------------//

inline ResponseIndex
ResponseIndexer::response_index(const size_t node_g,
                                const size_t surface,
                                const size_t index_s) const
{
  Require(node_g  < d_indices.size());
  Require(surface < d_indices[node_g].size());
  Require(index_s < d_indices[node_g][surface].size());
  return d_indices[node_g][surface][index_s];
}

inline ResponseIndex
ResponseIndexer::response_index(const size_t index_l) const
{
  Require(index_l < d_local_size);
  size_t node    = d_local_indices[index_l][0];
  size_t surface = d_local_indices[index_l][1];
  size_t nindex  = d_local_indices[index_l][2];
  return d_indices[node][surface][nindex];
}

//---------------------------------------------------------------------------//

inline ResponseIndexer::size_t ResponseIndexer::
nodal_to_local(const size_t node_g, const size_t index_n) const
{
  Require(node_g < d_number_local_nodes);
  return index_n + d_offsets[d_nodes->local_index(node_g)];
}

inline int ResponseIndexer::
global_to_local(const size_t index_g) const
{
  Require(index_g < d_global_size);
  int index_l = index_g - d_global_offset;
  if (index_g >= d_global_offset + d_local_size or
      index_g < d_global_offset)
  {
    index_l = -1;
  }
  return index_l;
}

inline ResponseIndexer::size_t ResponseIndexer::
nodal_to_global(const size_t node_g, const size_t index_n) const
{
  Require(node_g < d_indices.size());
  Require(index_n < d_sizes[node_g]);


  size_t index_g = index_n +
                   d_global_offsets[node_g];

  Ensure(index_g < d_global_size);
  return index_g;
}

inline ResponseIndexer::size_t ResponseIndexer::
local_to_global(const size_t index_l) const
{
  Require(index_l < d_local_size);
  size_t index_g = index_l + d_global_offset;
  Ensure(index_g < d_global_size);
  return index_g;
}

} // end namespace erme_response

#endif // RESPONSEINDEXER_I_HH_ 

//---------------------------------------------------------------------------//
//              end of file ResponseIndexer.i.hh
//---------------------------------------------------------------------------//
