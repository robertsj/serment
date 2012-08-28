//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   NodePartitioner.hh
 * \brief  NodePartitioner class definition
 * \author Jeremy Roberts
 * \date   Aug 24, 2012
 */
//---------------------------------------------------------------------------//

#ifndef NODEPARTITIONER_HH_
#define NODEPARTITIONER_HH_

#include "comm/Comm.hh"
#include "NodeList.hh"
#include <string>

namespace erme_geometry
{

/*!
 *  \class NodePartitioner
 *  \brief Partition nodes one level 1 communicator
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

  typedef serment_comm::Comm Comm;
  typedef Node::SP_node      SP_node;

  /*!
   *  \brief Constructor
   *  \param nodes  List of problem nodes
   */
  NodePartitioner();

  /// Partition the nodes
  void partition(NodeList &nodes);

private:

  /// Buffer for sending nodes
  std::string d_buffer;

  /// Buffer size
  int d_buffer_size;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /// Send a node from master to another process
  void broadcast_nodes(NodeList &nodes);

};

} // end namespace erme_geometry

#endif // NODEPARTITIONER_HH_ 

//---------------------------------------------------------------------------//
//              end of file NodePartitioner.hh
//---------------------------------------------------------------------------//
