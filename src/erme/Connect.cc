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
Connect::Connect(erme_geometry::NodeList &nodes,
                 erme_response::ResponseIndexer &indexer)
  : linear_algebra::Matrix(indexer.number_local_moments(),
                           indexer.number_local_moments(),
                           vec_int(indexer.number_local_moments(), 1),
                           vec_int(indexer.number_local_moments(), 1))
{
  using std::cout;
  using std::endl;

  Require(nodes.number_global_nodes() == indexer.number_nodes());

  // All nodes
  for (int n = nodes.lower_bound(); n < nodes.upper_bound(); n++)
  {
    int n_index = 0;

    // All surfaces
    for (int s = 0; s < nodes.node(n)->number_surfaces(); s++)
    {

      // Get the neighbor node and surface
      int neigh_n = nodes.neighbor(n, s).neighbor();
      int neigh_s = nodes.neighbor(n, s).surface();

      // Nodes should have consistent expansions on shared surfaces.  This
      // is a light check, since order counts can be the same with different
      // expansions.  A geometry preprocessor would be useful.
      if (neigh_n >= 0)
      {
        Assert(indexer.number_surface_moments(neigh_n, neigh_s)
               == indexer.number_surface_moments(n, s));
      }

      // All moments on a surface
      for (int m = 0; m < indexer.number_surface_moments(n, s); m++, n_index++)
      {
        // Set the row, column, and value
        int row      = indexer.global_index(n, n_index);
        int column   = indexer.global_index(n, n_index);
        double value = 0.0;

        if (nodes.neighbor(n, s).neighbor() >= 0)
        {
          // Get neighbor node moment index
          int neigh_n_index = indexer.node_index(neigh_n, neigh_s, m).local;

          // My row connects to their column.
          column = indexer.global_index(neigh_n, neigh_n_index);
          value  = 1.0;
        }
        else if (nodes.neighbor(n, s).neighbor() == erme_geometry::Node::REFLECT)
        {
          indexer.node_index(n, s, m).even_odd == 0 ? value = 1.0 : value = -1.0;
        }
        else if (nodes.neighbor(n, s).neighbor() == erme_geometry::Node::VACUUM)
        {
          // Nothing for vacuum
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

  // Assemble the matrix.
  assemble();
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file Connect.cc
//---------------------------------------------------------------------------//
