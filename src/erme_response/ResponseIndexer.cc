//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ResponseIndexer.cc
 *  @brief  ResponseIndexer
 *  @author Jeremy Roberts
 *  @date   Aug 24, 2012
 */
//---------------------------------------------------------------------------//

#include "ResponseIndexer.hh"
#include "comm/Comm.hh"
#include <iostream>

namespace erme_response
{

//---------------------------------------------------------------------------//
ResponseIndexer::ResponseIndexer(SP_db db, SP_nodelist nodes)
  : d_nodes(nodes)
  , d_order_reduction(0)
  , d_global_offset(0)
  , d_unique_size(0)
  , d_local_size(0)
  , d_global_size(0)
{
  // Preconditions
  Require(db);
  Require(d_nodes);
  // Dimension required for correct builder
  Insist(db->check("dimension"),
         "Parameter database must specify dimension.");
  int dimension = db->get<int>("dimension");
  Require(dimension >= 1);
  Require(dimension <= 3);
  // Check that all nodes have the same dimension
  for (int i = 0; i < d_nodes->number_unique_global_nodes(); ++i)
    Require(d_nodes->unique_node(i)->dimension() == dimension);

  // Initialize containers
  d_sizes.resize(d_nodes->number_unique_global_nodes(), 0);
  d_offsets.resize(d_nodes->number_local_nodes(), 0);
  d_global_offsets.resize(d_nodes->number_global_nodes(), 0);
  d_number_global_nodes = d_nodes->number_global_nodes();
  d_number_local_nodes = d_nodes->number_unique_local_nodes();

  // Get options from database.
  if (db->check("erme_order_reduction"))
    d_order_reduction = db->get<int>("erme_order_reduction");

  // Loop over all unique nodes.  All processes get all the same
  // info.  This makes it easier to ensure consistent expansions
  // between nodes.  The build routines return the size of the
  // unique node, which is recorded.
  for (size_t n = 0; n < d_nodes->number_unique_global_nodes(); n++)
  {
    // Moments size for the node.
    size_t size = 0;
    if (dimension == 1)
      size = build_1D(d_nodes->node(n), n);
    else if (dimension == 2)
      size = build_2D(d_nodes->node(n), n);
    else
      size = build_3D(d_nodes->node(n), n);
    d_sizes[n] = size;
    std::cout << " size of node " << n << " is " << size << std::endl;
  }

  /*
   *  Now d_sizes has the number of moments for each unique
   *  node.  Also, the indices for all the unique nodes have
   *  been built.  What we need now is to define the local
   *  and global numbers of moments and so forth that will be
   *  used in defining and updating global operators and
   *  vectors.
   */

  // Number of unique moments, i.e. the number we solve on local comm
  for (size_t node_ul = 0;
       node_ul < d_nodes->number_unique_local_nodes();
       ++node_ul)
  {
    // Get the unique global index.
    size_t node_ug  = d_nodes->unique_global_index_from_unique_local(node_ul);
    Assert(node_ug < d_sizes.size());
    // Add the size to the total unique number of moments.  This
    // is the number used to work share.
    d_unique_size += d_sizes[node_ug];
  }

  // Global size and offsets for each node
  for (size_t node_g = 0; node_g < d_nodes->number_global_nodes() - 1; ++node_g)
  {
    size_t node_ug = d_nodes->unique_global_index_from_global(node_g);
    d_global_offsets[node_g + 1] = d_global_offsets[node_g] + d_sizes[node_ug];
  }
  for (size_t node_g = 0; node_g < d_nodes->number_global_nodes(); ++node_g)
  {
    size_t node_ug = d_nodes->unique_global_index_from_global(node_g);
    d_global_size += d_sizes[node_ug];
  }

  // Compute local sizes, local offsets, and my global offset
  for (size_t node_g = nodes->lower_bound(); node_g < nodes->upper_bound(); ++node_g)
    d_local_size += d_sizes[d_nodes->unique_global_index_from_global(node_g)];
  for (int node_g = nodes->lower_bound(); node_g < nodes->upper_bound() - 1; ++node_g)
  {
    d_offsets[node_g - d_nodes->lower_bound() + 1] =
      d_offsets[node_g - d_nodes->lower_bound()] +
        d_sizes[d_nodes->unique_global_index_from_global(node_g)];
  }
  // My global offset (total number of moments before me)
  for (int node_g = 0; node_g < nodes->lower_bound(); ++node_g)
  {
    size_t node_ug = d_nodes->unique_global_index_from_global(node_g);
    d_global_offset += d_sizes[node_ug];
  }
  std::cout << " GLOBAL OFFSET = " << d_global_offset << std::endl;

  // Compute the unique local moment to (global node, surface, moment) index
  d_unique_indices.resize(d_unique_size, vec_size_t(3, 0));
  size_t unique_index = 0;
  for (int node_ul = 0; node_ul < nodes->number_unique_local_nodes(); ++node_ul)
  {
    size_t node_ug  = d_nodes->unique_global_index_from_unique_local(node_ul);
    for (int s = 0; s < d_indices[node_ug].size(); s++)
    {
      for (int m = 0; m < d_indices[node_ug][s].size(); m++, unique_index++)
      {
         d_unique_indices[unique_index][0] = node_ug;
         d_unique_indices[unique_index][1] = s;
         d_unique_indices[unique_index][2] = m;
      }
    }
  }

}

//---------------------------------------------------------------------------//
// Simple implementation for now.  A better format would be nice.
void ResponseIndexer::display() const
{
  using std::cout;
  using std::endl;
  cout << endl << "RESPONSE INDICES" << endl << endl;

  for (int node_g = 0; node_g < number_nodes(); ++node_g)
  {
    cout << "  global node " << node_g << endl;
    size_t node_ug = d_nodes->unique_global_index_from_global(node_g);
    Assert(node_ug < d_indices.size());
    for (int s = 0; s < d_indices[node_ug].size(); s++)
    {

      for (int m = 0; m < d_indices[node_ug][s].size(); m++)
      {
        cout << "    "
             << d_nodes->local_index_from_global(node_g)           << " "
             << response_index(node_ug, s, m).nodal    << " | "
             << response_index(node_ug, s, m).node     << " "
             << response_index(node_ug, s, m).surface  << " | "
             << response_index(node_ug, s, m).energy   << " | "
             << response_index(node_ug, s, m).polar    << " "
             << response_index(node_ug, s, m).azimuth  << " | "
             << response_index(node_ug, s, m).space0   << " "
             << response_index(node_ug, s, m).space1   << " | "
             << response_index(node_ug, s, m).even_odd << " |" << endl;
      } // end surface moment
    } // end surface
  } // end node
  cout << endl;
}


//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
ResponseIndexer::size_t
ResponseIndexer::build_3D(SP_node node, const size_t n)
{
  size_t nodal_index = 0;

  vec2_index surface_indices;
  for (size_t s = 0; s < node->number_surfaces(); s++)
  {
    // Max space order
    size_t mso = node->spatial_order(s, 0);
    if (node->spatial_order(s, 1) > mso) mso = node->spatial_order(s, 1);

    // Max angle order
    size_t mao = node->azimuthal_order(s);
    if (node->polar_order(s) > mao) mao = node->polar_order(s);

    // Max space-angle order
    size_t msao = mao;
    if (mso > mao) msao = mso;

    // Is this a south or north surface?
    bool south_north = (node->number_surfaces() - s > 2) ? false : true;

    vec_index moment_indices;
    for (size_t e = 0; e <= node->energy_order(s); e++)
    {
      for (size_t p = 0; p <= node->polar_order(s); p++)
      {
        for (size_t a = 0; a <= node->azimuthal_order(s); a++)
        {

          // Break if angle order too high
          size_t ao = p + a;
          if (d_order_reduction == 2 and ao > mao) break;

          for (size_t s0 = 0; s0 <= node->spatial_order(s, 0); s0++)
          {
            for (size_t s1 = 0; s1 <= node->spatial_order(s, 1); s1++)
            {

              // Break if space or space-angle order is too high
              size_t so = s0 + s1;
              if (d_order_reduction == 1 and so > mso) break;
              if (d_order_reduction == 3 and so + ao > msao) break;

              // Determine polarity. The polar polarity only contributes when
              // reflection is off a horizontal surface.  The other variables
              // always switch polarity.
              bool polarity =
                south_north ? (a + s0 + s1) % 2 : (a + p + s0 + s1) % 2;

              // Add the index
              moment_indices.push_back(
                ResponseIndex(n, s, e, p, a, s0, s1, polarity, nodal_index));

              // Update the local index
              nodal_index++;

            } // end second spatial
          } // end first spatial
        } // end azimuth
      } // end polar
    } // end energy
    surface_indices.push_back(moment_indices);
  } // end surface

  d_indices.push_back(surface_indices);
  return nodal_index;
}

//---------------------------------------------------------------------------//
ResponseIndexer::size_t
ResponseIndexer::build_2D(SP_node node, const size_t n)
{
  using std::cout;
  using std::endl;

  bool db = false;

  // Moment index local to a node
  size_t nodal_index = 0;

  if (db) std::cout << "NODE = " << n << std::endl;

  vec2_index surface_indices;
  for (size_t s = 0; s < node->number_surfaces(); s++)
  {
    if (db) cout << "SURFACE " << s << endl;

    // Max angle order
    size_t mao = node->azimuthal_order(s);
    if (node->polar_order(s) > mao) mao = node->polar_order(s);

    // Max space-angle order
    size_t msao = node->spatial_order(s, 0);
    if (mao > msao) msao = mao;

    vec_index moment_indices;
    for (size_t e = 0; e <= node->energy_order(s); e++)
    {
      if (db) cout << "  energy order " << e << endl;
      for (size_t p = 0; p <= node->polar_order(s); p++)
      {
        if (db) cout << "    polar order " << p << endl;
        for (size_t a = 0; a <= node->azimuthal_order(s); a++)
        {
          if (db) cout << "      azimuth order " << a << endl;
          // Break if angle order too high
          size_t ao = p + a;
          if (d_order_reduction == 2 and ao > mao) break;

          for (size_t s0 = 0; s0 <= node->spatial_order(s, 0); s0++)
          {
            if (db) cout << "        space order " << s0 << endl;
            // Break if space or space-angle order is too high
            if (d_order_reduction == 3 and s0 + ao > msao) break;

            // Determine polarity. The polar polarity only contributes when
            // reflection is off a horizontal surface.  The other variables
            // always switch polarity.
            bool polarity = (a + s0) % 2;
            if (db)  cout << "          polarity " << polarity << endl;
            // Add the index
            moment_indices.push_back(
              ResponseIndex(n, s, e, p, a, s0, 0, polarity, nodal_index));

            // Update the nodal index
            nodal_index++;

          } // end first spatial
        } // end azimuth
      } // end polar
    } // end energy
    surface_indices.push_back(moment_indices);
  } // end surface
  d_indices.push_back(surface_indices);

  return nodal_index;
}

//---------------------------------------------------------------------------//
ResponseIndexer::size_t
ResponseIndexer::build_1D(SP_node node, const size_t n)
{
  bool db = false;
  size_t nodal_index = 0;

  if(db) std::cout << " node  = " << n << std::endl;
  vec2_index surface_indices;
  for (size_t s = 0; s < node->number_surfaces(); s++)
  {
    if(db) std::cout << " surface  = " << s << std::endl;
    vec_index moment_indices;
    for (size_t e = 0; e <= node->energy_order(s); e++)
    {
      if(db) std::cout << " energy  = " << e << std::endl;
      for (size_t p = 0; p <= node->polar_order(s); p++)
      {
        if(db) std::cout << " polar  = " << p << std::endl;

        // Add the index
        moment_indices.push_back(
          ResponseIndex(n, s, e, p, 0, 0, 0, false, nodal_index));

        // Update the local index
        nodal_index++;

      } // end polar
    } // end energy
    surface_indices.push_back(moment_indices);
  } // end surface

  d_indices.push_back(surface_indices);

  return nodal_index;
}

} // end namespace erme_response

//---------------------------------------------------------------------------//
//              end of file ResponseIndexer.cc
//---------------------------------------------------------------------------//
