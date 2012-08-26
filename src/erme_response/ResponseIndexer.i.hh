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

inline ResponseIndexer::size_type
ResponseIndexer::number_nodes() const
{
  return d_sizes.size();
}

inline ResponseIndexer::size_type
ResponseIndexer::number_moments() const
{
  return d_local_size;
}

inline ResponseIndexer::size_type
ResponseIndexer::number_moments(const size_type node) const
{
  Require(node < d_sizes.size());
  return d_sizes[node];
}

inline ResponseIndexer::size_type
ResponseIndexer::global_number_moments() const
{
  return d_global_size;
}

inline ResponseIndex
ResponseIndexer::index(const size_type cindex, const size_type node) const
{
  Require(node < d_sizes.size());
  Require(cindex < d_sizes[node]);
  return d_indices[cindex + d_offsets[node]];
}

} // end namespace erme_response

#endif // RESPONSEINDEXER_I_HH_ 

//---------------------------------------------------------------------------//
//              end of file ResponseIndexer.i.hh
//---------------------------------------------------------------------------//
