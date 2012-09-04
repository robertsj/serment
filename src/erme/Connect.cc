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

/*
 *  Note, connectivity matrices have at most one nonzero value
 *  per row for normal use cases.  In the event that rotational
 *  symmetry is used, rows might have an additional entry, but
 *  the cost should not be great.
 *
 *  Also, since the connectivity matrix only really lives on the
 *  global communicator, the remaining processes end up doing
 *  wasted work.  Unless there is something better they can do,
 *  this should be fine.
 */
Connect::Connect(SP_nodelist nodes, SP_indexer indexer)
  : linear_algebra::Matrix(indexer->number_local_moments(),
                           indexer->number_local_moments(),
                           vec_int(indexer->number_local_moments(), 1),
                           vec_int(indexer->number_local_moments(), 1))
{
  using std::cout;
  using std::endl;

  Require(nodes->number_global_nodes() == indexer->number_nodes());

  // Loop over local nodes
  for (int n = nodes->lower_bound(); n < nodes->upper_bound(); n++)
  {

    int n_index = 0;

    // Loop over node surfaces
    for (int s = 0; s < nodes->node(n)->number_surfaces(); s++)
    {

      // Get the neighbor node and its surface
      int neigh_n = nodes->neighbor(n, s).neighbor();
      int neigh_s = nodes->neighbor(n, s).surface();

      // Nodes should have consistent expansions on shared surfaces.  This
      // is a light check, since order counts can be the same with different
      // expansions.  A geometry preprocessor would be useful.
      if (neigh_n >= 0)
      {
        Assert(indexer->number_surface_moments(neigh_n, neigh_s)
               == indexer->number_surface_moments(n, s));
      }

      // Loop over all moments on a surface
      for (int m = 0; m < indexer->number_surface_moments(n, s); m++, n_index++)
      {

        // Set the row, column, and value.  Note, if the surface is shared,
        // the column is changed below.
        int row      = indexer->nodal_to_global(n, n_index);
        int column   = indexer->nodal_to_global(n, n_index);
        double value = 0.0;

        if (nodes->neighbor(n, s).neighbor() >= 0)
        {
          // Get neighbor node moment index
          int neigh_n_index = indexer->response_index(neigh_n, neigh_s, m).local;

          // My row connects to their column.  Note, the matrix is symmetric
          // for most problems, but we let the other process fill the
          // reflection since it owns the row in which the entry lives.
          column = indexer->nodal_to_global(neigh_n, neigh_n_index);

//          if (serment_comm::Comm::rank() == 1)
//          {
//            cout << " neigh = " << neigh_n << " s = "
//                 << neigh_s << " m = " << m
//                 << " nidx = " << neigh_n_index
//                 << " col =  " << column
//                 << endl;
//          }

          value  = 1.0;
        }
        else if (nodes->neighbor(n, s).neighbor() == erme_geometry::Node::REFLECT)
        {
          // Reflection gets 1.0 or -1.0, depending on polarity of the moment
          indexer->response_index(n, s, m).even_odd
            == 0 ? value = 1.0 : value = -1.0;
        }
        else if (nodes->neighbor(n, s).neighbor() == erme_geometry::Node::VACUUM)
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
