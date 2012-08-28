//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseIndexer.hh
 * \brief  ResponseIndexer class definition
 * \author Jeremy Roberts
 * \date   Aug 24, 2012
 */
//---------------------------------------------------------------------------//

#ifndef RESPONSEINDEXER_HH_
#define RESPONSEINDEXER_HH_

#include "InputDB.hh"
#include "erme_geometry/NodeList.hh"
#include "ResponseIndex.hh"
#include <vector>

namespace erme_response
{

/*!
 *  \class ResponseIndexer
 *  \brief Indexes a node response vector
 *
 *  Each Node has assigned maximum orders in each phase space
 *  variable on each surface of the node.  Each variable is
 *  expanded in a set of basis functions, and so the combined
 *  response is a tensor product of such functions.  Hence,
 *  the total order of a particular response might be larger
 *  than the maximum allowed for any particular variable.  The
 *  user can limit the order of cross terms by setting the
 *  appropriate parameters in the database.
 *
 *  Note, the indexer is constructed by all processes and
 *  applies to the entire node list.
 *
 *  Relevant database entries:
 *    - dimension
 *    - erme_order_reduction
 */
class ResponseIndexer
{

public:

  typedef detran::InputDB::SP_input SP_db;
  typedef erme_geometry::NodeList::SP_node SP_node;
  typedef unsigned int size_type;
  typedef std::vector<ResponseIndex> vec_index;
  typedef std::vector<vec_index>     vec2_index;
  typedef std::vector<vec2_index>    vec3_index;
  /*!
   *  \brief Constructor
   *  \param node   Node to index
   */
  ResponseIndexer(SP_db db, erme_geometry::NodeList &nodes);

  /// Number of nodes indexed
  size_type number_nodes() const;

  /// Return the number of moments of a node
  size_type number_node_moments(const size_type node) const;

  /// Return the number of moments of a node surface
  size_type number_surface_moments(const size_type node,
                                   const size_type surface) const;

  /// Return the number of moments of all local nodes
  size_type number_local_moments() const;

  /// Return the number of moments of all nodes
  size_type number_global_moments() const;

  /*!
   *  \brief Get moment indices from cardinal index within node
   *  \param node     Global index of node
   *  \param surface  Surface index of node
   *  \param index    Moment index on surface of node
   */
  ResponseIndex node_index(const size_type node,
                           const size_type surface,
                           const size_type index) const;

  /// Get local moment index from a cardinal index within node
  size_type local_index(const size_type node, const size_type index) const
  {
    return index + d_offsets[node];
  }

  /// Get global moment index from a cardinal index within node
  size_type global_index(const size_type node, const size_type index) const
  {
    return index + d_offsets[node] + d_global_offset;
  }

  /// Display the indices in a nice format
  void display() const;

private:

  //---------------------------------------------------------------------------//
  // PRIVATE DATA
  //---------------------------------------------------------------------------//

  /// List of indices for all nodes [nodes][surface][moments]
  vec3_index d_indices;

  /// Vector of moment sizes per node
  std::vector<size_type> d_sizes;

  /// Offset of node indices (within local set)
  std::vector<size_type> d_offsets;

  /// Local size of moments vector
  size_type d_local_size;

  /// Global offset (i.e. the number of moments before me)
  size_type d_global_offset;

  /// Global size
  size_type d_global_size;

  /*!
   *  \brief Order reduction selector
   *
   *  In general, each phase space variable is treated independently via
   *  a one dimensional orthogonal basis set.  The full phase space is
   *  approximated via a tensor product of these one dimensional functions.
   *
   *  To limit the order of cross terms, the user can select an order
   *  reduction approach.  The options currently implemented are
   *
   *    0 - no reduction
   *    1 - spatial order is limited by the largest of the two spatial
   *        orders (applicable to 3D only)
   *    2 - angular order is limited by the largest of azimuthal and
   *        polar order (not applicable to 1D)
   *    3 - combination of 1 and 2
   *    4 - maximum space-angle order is the maximum of the available
   *        space and angle orders (not applicable to 1D)
   *
   *  Note, this reduction scheme is applied on each surface
   *  seperately.
   */
  int d_order_reduction;

  //---------------------------------------------------------------------------//
  // IMPLEMENTATION
  //---------------------------------------------------------------------------//

  // Build the indices for a node and return the number of moments
  size_type build_1D(SP_node node, const size_type n);
  size_type build_2D(SP_node node, const size_type n);
  size_type build_3D(SP_node node, const size_type n);

};

} // end namespace erme_response

// Inline member definitions
#include "ResponseIndexer.i.hh"

#endif // RESPONSEINDEXER_HH_ 

//---------------------------------------------------------------------------//
//              end of file ResponseIndexer.hh
//---------------------------------------------------------------------------//
