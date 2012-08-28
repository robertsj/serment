//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Connect.cc
 * \brief  Connect 
 * \author Jeremy Roberts
 * \date   Aug 23, 2012
 */
//---------------------------------------------------------------------------//

#include "Connect.hh"
#include <iostream>

namespace erme
{

// Note, connectivity matrices have at most one nonzero value
// per row.
Connect::Connect(erme_geometry::NodeList &list,
                 erme_response::ResponseIndexer &indexer)
  : linear_algebra::Matrix(indexer.number_local_moments(),
                           indexer.number_local_moments(),
                           vec_int(indexer.number_local_moments(), 1),
                           vec_int(indexer.number_local_moments(), 1))
{
  using std::cout;
  using std::endl;

  Require(list.number_global_nodes() == indexer.number_nodes());

  // All nodes
  for (int n = list.lower_bound(); n < list.upper_bound(); n++)
  {
    cout << "n = " << n << endl;
    int idx = 0;
    // All surfaces
    for (int s = 0; s < list.node(n)->number_surfaces(); s++)
    {
      cout << "  s = " << s << endl;

      // Get the neighbor node and surface
      int neigh_n = list.neighbor(n, s).neighbor();
      int neigh_s = list.neighbor(n, s).surface();


      // Nodes should have consistent expansions on shared surfaces
      if (neigh_n >= 0)
      {
        Assert(indexer.number_surface_moments(neigh_n, neigh_s)
               == indexer.number_surface_moments(n, s));
      }
      // All moments
      for (int m = 0; m < indexer.number_surface_moments(n, s); m++, idx++)
      {
        cout << "    m = " << m << endl;

        // Set the row, column, and value
        int row    = indexer.global_index(n, idx);
        int column = indexer.global_index(n, idx);

        cout << "      row = " << row << endl;
        double value = 0.0;
        if (list.neighbor(n, s).neighbor() >= 0)
        {
          // Connect my row to their column
          row    = indexer.global_index(n, idx);
          column = indexer.global_index(n, idx);
          value  = 2.0;
        }
        else if (list.neighbor(n, s).neighbor() == erme_geometry::Node::REFLECT)
        {
          indexer.node_index(n, s, m).even_odd == 0 ? value = 1.0 : value = -1.0;
        }
        else if (list.neighbor(n, s).neighbor() == erme_geometry::Node::VACUUM)
        {
          // vacuum gets zero
        }
        else
        {
          THROW("Invalid neighbor value");
        }

        // Insert the value.
        insert_values(1, &row, 1, &column, &value);

      } // end moments
    } // end surfaces
  } // end nodes


}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file Connect.cc
//---------------------------------------------------------------------------//
