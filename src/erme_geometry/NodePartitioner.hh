//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   NodePartitioner.hh
 *  @brief  NodePartitioner class definition
 *  @author Jeremy Roberts
 *  @date   Aug 24, 2012
 */
//---------------------------------------------------------------------------//

#ifndef erme_geometry_NODEPARTITIONER_HH_
#define erme_geometry_NODEPARTITIONER_HH_

#include "comm/Comm.hh"
#include "NodeList.hh"
#include <string>

namespace erme_geometry
{

/**
 *  @class NodePartitioner
 *  @brief Partition nodes one level 1 communicator
 *
 *  This is a very light weight partitioning that simply broadcasts
 *  the list of nodes and assigns array bounds for each receiving
 *  process.
 */
class NodePartitioner
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef serment_comm::Comm      Comm;
  typedef Node::SP_node           SP_node;
  typedef NodeList::SP_nodelist   SP_nodelist;
  typedef unsigned int            size_t;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /// Constructor
  NodePartitioner();

  /**
   *  @brief Partition a list of nodes
   *  @param nodes  Node list
   */
  void partition(SP_nodelist &nodes);

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Buffer for sending nodes
  std::string d_buffer;

  /// Buffer size
  int d_buffer_size;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /// Send a node from master to another process
  void broadcast_nodes(SP_nodelist &nodes);

};

} // end namespace erme_geometry

#endif // erme_geometry_NODEPARTITIONER_HH_

//---------------------------------------------------------------------------//
//              end of file NodePartitioner.hh
//---------------------------------------------------------------------------//
