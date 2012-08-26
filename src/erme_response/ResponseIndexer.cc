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

namespace erme_response
{

ResponseIndexer::ResponseIndexer(SP_db db, erme_geometry::NodeList &nodes)
  : d_indices(nodes.number_nodes())
  , d_sizes(nodes.number_nodes(), 0)
  , d_offsets(nodes.number_nodes(), 0)
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
  for (size_type n = 0; n < nodes.number_nodes(); n++)
  {
    size_type size = 0;
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

}

ResponseIndexer::size_type
ResponseIndexer::build_3D(SP_node node, const size_type n)
{

  size_type local_index = 0;
  size_type offset;
  for (size_type i = 0; i < n; i++)
    offset += d_sizes[i];

  for (size_type s = 0; s < node->number_surfaces(); s++)
  {
    // Max space order
    size_type mso = node->spatial_order(s, 0);
    if (node->spatial_order(s, 1) > mso) mso = node->spatial_order(s, 1);

    // Max angle order
    size_type mao = node->azimuthal_order(s);
    if (node->polar_order(s) > mao) mao = node->polar_order(s);

    // Max space-angle order
    size_type msao = mao;
    if (mso > mao) msao = mso;

    // Is this a south or north surface?
    bool south_north = (node->number_surfaces() - s > 2) ? false : true;

    for (size_type e = 0; e < node->energy_order(s); e++)
    {
      for (size_type p = 0; p < node->polar_order(s); p++)
      {
        for (size_type a = 0; a < node->azimuthal_order(s); a++)
        {

          // Break if angle order too high
          size_type ao = p + a;
          if (d_order_reduction == 2 and ao > mao) break;

          for (size_type s0 = 0; s0 < node->spatial_order(s, 0); s0++)
          {
            for (size_type s1 = 0; s1 < node->spatial_order(s, 1); s1++)
            {

              // Break if space or space-angle order is too high
              size_type so = s0 + s1;
              if (d_order_reduction == 1 and so > mso) break;
              if (d_order_reduction == 3 and so + ao > msao) break;

              // Determine polarity. The polar polarity only contributes when
              // reflection is off a horizontal surface.  The other variables
              // always switch polarity.
              bool polarity =
                south_north ? (a + s0 + s1) % 2 : (a + p + s0 + s1) % 2;

              // Add the index
              ResponseIndex tmp(s, e, p, a, s0, s1, polarity);
              d_indices[local_index + offset] = tmp;
              // ResponseIndex(s, e, p, a, s0, 0, polarity);

              // Update the local index
              local_index++;

            } // end second spatial
          } // end first spatial
        } // end azimuth
      } // end polar
    } // end energy
  } // end surface

  return local_index;
}

ResponseIndexer::size_type
ResponseIndexer::build_2D(SP_node node, const size_type n)
{

  size_type local_index = 0;
  size_type offset;
  for (size_type i = 0; i < n; i++)
    offset += d_sizes[i];

  for (size_type s = 0; s < node->number_surfaces(); s++)
  {

    // Max angle order
    size_type mao = node->azimuthal_order(s);
    if (node->polar_order(s) > mao) mao = node->polar_order(s);

    // Max space-angle order
    size_type msao = node->spatial_order(s, 0);
    if (msao > mao) msao = mao;

    for (size_type e = 0; e < node->energy_order(s); e++)
    {
      for (size_type p = 0; p < node->polar_order(s); p++)
      {
        for (size_type a = 0; a < node->azimuthal_order(s); a++)
        {

          // Break if angle order too high
          size_type ao = p + a;
          if (d_order_reduction == 2 and ao > mao) break;

          for (size_type s0 = 0; s0 < node->spatial_order(s, 0); s0++)
          {

            // Break if space or space-angle order is too high
            if (d_order_reduction == 3 and s0 + ao > msao) break;

            // Determine polarity. The polar polarity only contributes when
            // reflection is off a horizontal surface.  The other variables
            // always switch polarity.
            bool polarity = (a + s0) % 2;

            // Add the index
            ResponseIndex tmp(s, e, p, a, s0, 0, polarity);
            d_indices[local_index + offset] = tmp;
            // ResponseIndex(s, e, p, a, s0, 0, polarity);

            // Update the local index
            local_index++;

          } // end first spatial
        } // end azimuth
      } // end polar
    } // end energy
  } // end surface

  return local_index;
}

ResponseIndexer::size_type
ResponseIndexer::build_1D(SP_node node, const size_type n)
{

  size_type local_index = 0;
  size_type offset;
  for (size_type i = 0; i < n; i++)
    offset += d_sizes[i];

  for (size_type s = 0; s < node->number_surfaces(); s++)
  {
    for (size_type e = 0; e < node->energy_order(s); e++)
    {
      for (size_type p = 0; p < node->polar_order(s); p++)
      {

        // Add the index
        ResponseIndex tmp(s, e, p, 0, 0, 0, false);
        d_indices[local_index + offset] = tmp;
        // ResponseIndex(s, e, p, a, s0, 0, polarity);

        // Update the local index
        local_index++;

      } // end polar
    } // end energy
  } // end surface

  return local_index;
}

} // end namespace erme_response

//---------------------------------------------------------------------------//
//              end of file ResponseIndexer.cc
//---------------------------------------------------------------------------//
