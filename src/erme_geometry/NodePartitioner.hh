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
// Serialization to string
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/device/array.hpp>
// Serialization
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

namespace erme_geometry
{

/*!
 *  \class NodePartitioner
 *  \brief Partition nodes across level 1 process groups
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
  void send_node(SP_node node, int destination);

  /// Receive a node from master
  void receive_node(SP_node &node);

};

} // end namespace erme_geometry

#endif // NODEPARTITIONER_HH_ 

//---------------------------------------------------------------------------//
//              end of file NodePartitioner.hh
//---------------------------------------------------------------------------//
