//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseIndexer.cc
 * \brief  ResponseIndexer 
 * \author Jeremy Roberts
 * \date   Aug 24, 2012
 */
//---------------------------------------------------------------------------//

#include "ResponseIndexer.hh"
#include "comm/Comm.hh"
#include <iostream>
namespace erme_response
{

ResponseIndexer::ResponseIndexer(SP_db db, erme_geometry::NodeList &nodes)
  : d_sizes(nodes.number_global_nodes(), 0)
  , d_offsets(nodes.number_global_nodes(), 0)
  , d_order_reduction(0)
  , d_local_size(0)
  , d_global_offset(0)
  , d_global_size(0)
{
  // Preconditions
  Require(db);

  Insist(db->check("dimension"), "Parameter database must specify dimension.");
  int dimension = db->get<int>("dimension");

  // Get options from database.
  if (db->check("erme_order_reduction"))
    d_order_reduction = db->get<int>("erme_order_reduction");

  // Loop over all nodes
  for (size_t n = 0; n < nodes.number_global_nodes(); n++)
  {
    // Moments size for the node
    size_t size = 0;
    if (dimension == 1)
      size = build_1D(nodes.node(n), n);
    else if (dimension == 2)
      size = build_2D(nodes.node(n), n);
    else
      size = build_3D(nodes.node(n), n);
    d_sizes[n] = size;
    d_offsets[n] = d_local_size;
    d_local_size += size;
  }

  // Communicate local sizes
  std::vector<int> all_sizes(serment_comm::Comm::size(), 0);
  all_sizes[serment_comm::Comm::rank()] = d_local_size;
  serment_comm::Comm::global_sum(&all_sizes[0], all_sizes.size());

  // Compute the global moment offset
  for (int i = 0; i < serment_comm::Comm::rank(); i++)
    d_global_offset += all_sizes[i];

  // Compute the global size
  for (int i = 0; i < all_sizes.size(); i++)
    d_global_size += all_sizes[i];

  // Compute the local moment to (node, surface, moment) index
  d_local_indices.resize(d_local_size, vec_size_t(3, 0));
  size_t local_index = 0;
  for (int n = nodes.lower_bound(); n < nodes.upper_bound(); n++)
  {
    for (int s = 0; s < d_indices[n].size(); s++)
    {
      for (int m = 0; m < d_indices[n][s].size(); m++, local_index++)
      {
         d_local_indices[local_index][0] = n;
         d_local_indices[local_index][0] = s;
         d_local_indices[local_index][0] = m;
      }
    }
  }

}

// Simple implementation for now.  A better format would be nice.
void ResponseIndexer::display() const
{
  using std::cout;
  using std::endl;
  for (int n = 0; n < number_nodes(); n++)
  {
    cout << "GLOBAL NODE " << n << endl;
    for (int s = 0; s < d_indices[n].size(); s++)
    {
      for (int m = 0; m < d_indices[n][s].size(); m++)
      {
        cout << "    " << m << " | "
             << node_index(n, m, s).surface  << " | "
             << node_index(n, m, s).energy   << " | "
             << node_index(n, m, s).polar    << " "
             << node_index(n, m, s).azimuth  << " | "
             << node_index(n, m, s).space0   << " "
             << node_index(n, m, s).space1   << " | "
             << node_index(n, m, s).even_odd << " |" << endl;
      }
    }
  }
}


//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//

ResponseIndexer::size_t
ResponseIndexer::build_3D(SP_node node, const size_t n)
{

  size_t local_index = 0;
  size_t offset;
  for (size_t i = 0; i < n; i++)
    offset += d_sizes[i];

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
              ResponseIndex tmp(s, e, p, a, s0, s1, polarity, local_index);
              moment_indices.push_back(tmp);

              // Update the local index
              local_index++;

            } // end second spatial
          } // end first spatial
        } // end azimuth
      } // end polar
    } // end energy
    surface_indices.push_back(moment_indices);
  } // end surface

  d_indices.push_back(surface_indices);
  return local_index;
}

ResponseIndexer::size_t
ResponseIndexer::build_2D(SP_node node, const size_t n)
{
  using std::cout;
  using std::endl;

  // Moment index local to a node
  size_t local_index = 0;
  // Number of moments on preceding nodes
  size_t offset = 0;
  for (size_t i = 0; i < n; i++)
    offset += d_sizes[i];

  vec2_index surface_indices;
  for (size_t s = 0; s < node->number_surfaces(); s++)
  {
    //cout << "SURFACE " << s << endl;

    // Max angle order
    size_t mao = node->azimuthal_order(s);
    if (node->polar_order(s) > mao) mao = node->polar_order(s);

    // Max space-angle order
    size_t msao = node->spatial_order(s, 0);
    if (mao > msao) msao = mao;

    vec_index moment_indices;
    for (size_t e = 0; e <= node->energy_order(s); e++)
    {
      //cout << "  energy order " << e << endl;
      for (size_t p = 0; p <= node->polar_order(s); p++)
      {
        //cout << "    polar order " << p << endl;
        for (size_t a = 0; a <= node->azimuthal_order(s); a++)
        {
          //cout << "      azimuth order " << a << endl;
          // Break if angle order too high
          size_t ao = p + a;
          if (d_order_reduction == 2 and ao > mao) break;

          for (size_t s0 = 0; s0 <= node->spatial_order(s, 0); s0++)
          {
            //cout << "        space order " << s0 << endl;
            // Break if space or space-angle order is too high
            if (d_order_reduction == 3 and s0 + ao > msao) break;

            // Determine polarity. The polar polarity only contributes when
            // reflection is off a horizontal surface.  The other variables
            // always switch polarity.
            bool polarity = (a + s0) % 2;

            // Add the index
            moment_indices.push_back(ResponseIndex(s, e, p, a, s0, 0, polarity, local_index));

            // Update the local index
            local_index++;

          } // end first spatial
        } // end azimuth
      } // end polar
    } // end energy
    surface_indices.push_back(moment_indices);
  } // end surface
  d_indices.push_back(surface_indices);

  return local_index;
}

ResponseIndexer::size_t
ResponseIndexer::build_1D(SP_node node, const size_t n)
{

  size_t local_index = 0;
  size_t offset;
  for (size_t i = 0; i < n; i++)
    offset += d_sizes[i];

  vec2_index surface_indices;
  for (size_t s = 0; s <= node->number_surfaces(); s++)
  {
    vec_index moment_indices;
    for (size_t e = 0; e <= node->energy_order(s); e++)
    {
      for (size_t p = 0; p <= node->polar_order(s); p++)
      {

        // Add the index
        moment_indices.push_back(ResponseIndex(s, e, p, 0, 0, 0, false, local_index));

        // Update the local index
        local_index++;

      } // end polar
    } // end energy
    surface_indices.push_back(moment_indices);
  } // end surface

  d_indices.push_back(surface_indices);

  return local_index;
}

} // end namespace erme_response

//---------------------------------------------------------------------------//
//              end of file ResponseIndexer.cc
//---------------------------------------------------------------------------//
