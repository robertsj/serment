//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ResponseIndexer.i.hh
 *  @brief  ResponseIndexer inline member definitions
 *  @author Jeremy Roberts
 *  @date   Aug 24, 2012
 */
//---------------------------------------------------------------------------//

#ifndef erme_response_RESPONSEINDEXER_I_HH_
#define erme_response_RESPONSEINDEXER_I_HH_

namespace erme_response
{

//---------------------------------------------------------------------------//
inline ResponseIndexer::size_t
ResponseIndexer::number_nodes() const
{
  return d_nodes->number_global_nodes();
}

//---------------------------------------------------------------------------//
inline ResponseIndexer::size_t
ResponseIndexer::number_surface_moments(const size_t node_ug,
                                        const size_t surface) const
{
  Require(node_ug < d_indices.size());
  Require(surface < d_indices[node_ug].size());
  return d_indices[node_ug][surface].size();
}

//---------------------------------------------------------------------------//
inline ResponseIndexer::size_t
ResponseIndexer::number_node_moments(const size_t node_ug) const
{
  // Preconditions
  Require(node_ug < d_sizes.size());

  return d_sizes[node_ug];
}

//---------------------------------------------------------------------------//
inline ResponseIndexer::size_t
ResponseIndexer::number_unique_moments() const
{
  return d_unique_size;
}
//---------------------------------------------------------------------------//
inline ResponseIndexer::size_t
ResponseIndexer::number_local_moments() const
{
  return d_local_size;
}
//---------------------------------------------------------------------------//
inline ResponseIndexer::size_t
ResponseIndexer::number_global_moments() const
{
  return d_global_size;
}

//---------------------------------------------------------------------------//
inline ResponseIndex
ResponseIndexer::response_index(const size_t node_ug,
                                const size_t surface,
                                const size_t index_s) const
{
  Require(node_ug < d_indices.size());
  Require(surface < d_indices[node_ug].size());
  Require(index_s < d_indices[node_ug][surface].size());
  return d_indices[node_ug][surface][index_s];
}

//---------------------------------------------------------------------------//
inline ResponseIndex
ResponseIndexer::response_index_from_unique_local(const size_t index_ul) const
{
  Require(index_ul < d_unique_size);

  size_t node_ug = d_unique_indices[index_ul][0];
  size_t surface = d_unique_indices[index_ul][1];
  size_t nindex  = d_unique_indices[index_ul][2];

  Ensure(node_ug < d_indices.size());
  Ensure(surface < d_indices[node_ug].size());
  Ensure(nindex  < d_indices[node_ug][surface].size());
  return d_indices[node_ug][surface][nindex];
}

//---------------------------------------------------------------------------//
inline ResponseIndex
ResponseIndexer::response_index_from_local(const size_t index_l) const
{
  Require(index_l < number_local_moments());

  size_t index_ul = local_index_to_unique(index_l);

  Ensure(index_ul < d_unique_size);
  return response_index_from_unique_local(index_ul);
}

//---------------------------------------------------------------------------//
inline ResponseIndexer::size_t ResponseIndexer::
nodal_index_to_local(const size_t node_g, const size_t index_n) const
{
  Require(node_g < d_number_local_nodes);
  return index_n + d_offsets[d_nodes->local_index_from_global(node_g)];
}

//---------------------------------------------------------------------------//
inline int ResponseIndexer::
global_index_to_local(const size_t index_g) const
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

//---------------------------------------------------------------------------//
inline ResponseIndexer::size_t ResponseIndexer::
nodal_index_to_global(const size_t node_g, const size_t index_n) const
{
  Require(node_g < d_global_offsets.size());
  size_t node_ug = d_nodes->unique_global_index_from_global(node_g);
  Require(node_ug < d_indices.size());
  Require(index_n < d_sizes[node_ug]);

  size_t index_g = index_n + d_global_offsets[node_g];

  Ensure(index_g < d_global_size);
  return index_g;
}

//---------------------------------------------------------------------------//
inline ResponseIndexer::size_t ResponseIndexer::
local_index_to_global(const size_t index_l) const
{
  Require(index_l < d_local_size);

  size_t index_g = index_l + d_global_offset;

  Ensure(index_g < d_global_size);
  return index_g;
}

//---------------------------------------------------------------------------//
inline ResponseIndexer::size_t ResponseIndexer::
local_index_to_unique(const size_t index_l) const
{
  Require(index_l < d_local_to_unique.size());

  size_t index_ul = d_local_to_unique[index_l];

  Ensure(index_ul < number_unique_moments());
  return index_ul;
}

} // end namespace erme_response

#endif // erme_responseRESPONSEINDEXER_I_HH_

//---------------------------------------------------------------------------//
//              end of file ResponseIndexer.i.hh
//---------------------------------------------------------------------------//
