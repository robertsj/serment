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
Connect::Connect(erme_geometry::NodeList &list,
                 erme_response::ResponseIndexer &indexer)
  : linear_algebra::Matrix(indexer.number_moments(),
                           indexer.number_moments(),
                           vec_int(indexer.number_moments(), 1),
                           vec_int(indexer.number_moments(), 1))
{

  Require(list.number_nodes() == indexer.number_nodes());

  for (int n = 0; n < indexer.number_nodes(); n++)
  {
    for (int m = 0; m < indexer.number_moments(n); m++)
    {
      int s    = indexer.index(m, n).surface;
      int s_neigh = list.neighbor(n, s);

      if (list.neighbor(n, s) == erme_geometry::Node::VACUUM)
      {

      }


//      insert_values(const size_type number_rows,
//                    const int *rows,
//                    const size_type number_columns,
//                    const int *columns,
//                    const double *values);

    }
  }


}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file Connect.cc
//---------------------------------------------------------------------------//
