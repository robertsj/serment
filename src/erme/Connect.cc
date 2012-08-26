//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Connect.cc
 * \brief  Connect 
 * \author Jeremy Roberts
 * \date   Aug 23, 2012
 */
//---------------------------------------------------------------------------//

#include "Connect.hh"

namespace erme
{

// Note, connectivity matrices have at most one nonzero value
// per row.
Connect::Connect(const erme_geometry::NodeList &list,
                 const erme_response::ResponseIndexer &indexer)
  : linear_algebra::Matrix(indexer.number_moments(),
                           indexer.number_moments(),
                           vec_int(indexer.number_moments(), 1),
                           vec_int(indexer.number_moments(), 1))
{
  /* ... */
}


} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file Connect.cc
//---------------------------------------------------------------------------//
