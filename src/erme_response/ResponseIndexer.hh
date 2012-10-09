//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ResponseIndexer.hh
 *  @brief  ResponseIndexer class definition
 *  @author Jeremy Roberts
 *  @date   Aug 24, 2012
 */
//---------------------------------------------------------------------------//

#ifndef RESPONSEINDEXER_HH_
#define RESPONSEINDEXER_HH_

#include "ResponseIndex.hh"
#include "erme_geometry/NodeList.hh"
#include "utilities/InputDB.hh"
#include "utilities/SP.hh"
#include <vector>

namespace erme_response
{

/**
 *  @class ResponseIndexer
 *  @brief Indexes a node response vector
 *
 *  Each Node has assigned maximum orders in each phase space
 *  variable on each surface of the node.  Each variable is
 *  expanded in a set of basis functions, and so the combined
 *  response is a tensor product of such functions.  Hence,
 *  the total order of a particular response might be larger
 *  than the maximum allowed for any particular variable.  The
 *  user can limit the order of cross terms by setting the
 *  order reduction level.  The current options are:
 *    0 - no reduction (default)
 *    1 - spatial order is limited by the largest of the two spatial
 *        orders (applicable to 3D only)
 *    2 - angular order is limited by the largest of azimuthal and
 *        polar order (not applicable to 1D)
 *    3 - combination of 1 and 2
 *    4 - maximum space-angle order is the maximum of the available
 *        space and angle orders (not applicable to 1D)
 *
 *  Note, the indexer is constructed by all processes and
 *  applies to the entire node list.  However, a number of index
 *  functions help provide indices into the global, local, and
 *  nodal moments.
 *
 *  Relevant database entries:
 *    - dimension
 *    - erme_order_reduction
 */
/**
 *  @example erme_response/test/test_ResponseIndexer.cc
 *
 *  Test of ResponseIndexer class.
 */
class ResponseIndexer
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<ResponseIndexer> SP_indexer;
  typedef detran_utilities::InputDB::SP_input   SP_db;
  typedef erme_geometry::NodeList::SP_nodelist  SP_nodelist;
  typedef erme_geometry::NodeList::SP_node      SP_node;
  typedef detran_utilities::size_t              size_t;
  typedef detran_utilities::vec_size_t          vec_size_t;
  typedef detran_utilities::vec2_size_t         vec2_size_t;
  typedef std::vector<ResponseIndex>            vec_index;
  typedef std::vector<vec_index>                vec2_index;
  typedef std::vector<vec2_index>               vec3_index;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param db     Pointer to parameter database
   *  @param nodes  Pointer to node list
   */
  ResponseIndexer(SP_db db, SP_nodelist nodes);

  //-------------------------------------------------------------------------//

  /// Total number of nodes in the problem
  size_t number_nodes() const;

  /**
   *  @brief Number of moments associated with a node
   *  @param node_g   Global node index
   */
  size_t number_node_moments(const size_t node_g) const;

  /**
   *  @brief Number of moments on a node surface
   *  @param node_g   Global node index
   *  @param surface  Node surface index
   */
  size_t number_surface_moments(const size_t node_g,
                                const size_t surface_n) const;

  /// Return the number of moments of all local nodes
  size_t number_local_moments() const;

  /// Return the number of moments of all nodes
  size_t number_global_moments() const;

  //-------------------------------------------------------------------------//

  /**
   *  @brief Get moment indices from cardinal index within node
   *  @param node_g       Global index of node
   *  @param surface_n    Surface index of node
   *  @param index_s      Moment index on surface of node
   */
  ResponseIndex response_index(const size_t node_g,
                               const size_t surface_n,
                               const size_t index_s) const;

  /**
   *  @brief Get moment indices from local cardinal index
   *  @param index_l      Moment index within local nodes
   */
  ResponseIndex response_index(const size_t index_l) const;

  //-------------------------------------------------------------------------//

  /**
   *  @brief Get local moment index from a cardinal index within node
   *  @param  node_l      Local node index
   *  @param  index_n     Cardinal moment index within node
   */
  size_t nodal_to_local(const size_t node_l, const size_t index_n) const;

  /**
   *  @brief Get local moment index from global moment index
   *
   *  Note, this returns -1 if the global index doesn't represent
   *  a local index.
   *
   *  @param  index_g   Global node index
   */
  int global_to_local(const size_t index_g) const;

  /**
   *  @brief Get global moment index from a nodal moment index
   *  @param  node_g    Global node index
   *  @param  index_n   Moment index within node
   */
  size_t nodal_to_global(const size_t node_g, const size_t index_n) const;

  /**
   *  @brief Get global moment index from a local moment index
   *  @param  index_l   Local moment index
   */
  size_t local_to_global(const size_t index_l) const;

  //-------------------------------------------------------------------------//

  /// Display the indices in a nice format
  void display() const;

private:

  //-------------------------------------------------------------------------//
  // PRIVATE DATA
  //-------------------------------------------------------------------------//

  /// Node list
  SP_nodelist d_nodes;
  /// List of indices for all nodes [nodes][surface][moments]
  vec3_index d_indices;
  /// Map a local cardinal index to the [node, surface, moment]
  vec2_size_t d_local_indices;
  /// Vector of moment sizes per node
  std::vector<size_t> d_sizes;
  /// Offset of node indices (within local set)
  std::vector<size_t> d_offsets;
  /// Offset of node indices (within local set)
  std::vector<size_t> d_global_offsets;
  /// Local size of moments vector
  size_t d_local_size;
  /// Global offset (i.e. the number of moments before me)
  size_t d_global_offset;
  /// Global size
  size_t d_global_size;
  /// Number global nodes
  size_t d_number_global_nodes;
  /// Number local nodes
  size_t d_number_local_nodes;
  /// Order reduction selector
  int d_order_reduction;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  // Build the indices for a node and return the number of moments
  size_t build_1D(SP_node node, const size_t n);
  size_t build_2D(SP_node node, const size_t n);
  size_t build_3D(SP_node node, const size_t n);

};

} // end namespace erme_response

//---------------------------------------------------------------------------//
// INLINE MEMBER DEFINITIONS
//---------------------------------------------------------------------------//

#include "ResponseIndexer.i.hh"

#endif // RESPONSEINDEXER_HH_ 

//---------------------------------------------------------------------------//
//              end of file ResponseIndexer.hh
//---------------------------------------------------------------------------//
