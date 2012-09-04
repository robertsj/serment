//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   NodePartitioner.cc
 * \brief  NodePartitioner member definitions
 * \author Jeremy Roberts
 * \date   Aug 24, 2012
 */
//---------------------------------------------------------------------------//

#include "NodePartitioner.hh"

#include <iostream>
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

/*
 *  This is a simple partitioning process which merely
 *  broadcasts the node list to all processes, and assigns
 *  to each process the proper bounds of the global node
 *  array.
 *
 */

NodePartitioner::NodePartitioner()
  : d_buffer_size(0)
{
  /* ... */
}

void NodePartitioner::partition(SP_nodelist &nodes)
{

  size_t number_nodes = 0;
  size_t local_number_nodes = 0;
  size_t lower_bound = 0;

  // Switch to global communicator
  if (Comm::is_global())
  {
    Comm::set(serment_comm::global);
    // Partition the nodes using Comm's default partitioner
    if (Comm::rank() == 0)
    {
      number_nodes = nodes->number_global_nodes();
    }
    Comm::partition(number_nodes, lower_bound, local_number_nodes);
  }

  // Switch to local and broadcast local bound and size
  Comm::set(serment_comm::local);
  Comm::broadcast(&lower_bound, 1, 0);
  Comm::broadcast(&local_number_nodes, 1, 0);

  // Switch to world, broadcast nodes to everyone, and finalize.
  Comm::set(serment_comm::world);
  broadcast_nodes(nodes);
  nodes->set_bounds(lower_bound, lower_bound + local_number_nodes);
  nodes->finalize();

}

void NodePartitioner::broadcast_nodes(SP_nodelist &nodes)
{
  // Preconditions
  if (Comm::rank() == 0)
  {
    Require(nodes);
  }
  Require(serment_comm::communicator == serment_comm::world);

  // Clear the buffer.  This wipes the contents, but keeps the
  // memory allocated.
  d_buffer.clear();

  if (Comm::rank() == 0)
  {
    // Setup buffer and archive
    typedef boost::iostreams::back_insert_device<std::string> insert_t;
    insert_t inserter(d_buffer);
    boost::iostreams::stream<insert_t> s(inserter);
    boost::archive::binary_oarchive output_archive(s);

    // Archive the node and neighbor
    output_archive << nodes;
    s.flush();

    // Send the archive
    d_buffer_size = d_buffer.size();
    Comm::broadcast(&d_buffer_size, 1, 0);
    Comm::broadcast((char*)d_buffer.data(), d_buffer_size, 0);
  }
  else
  {
    // Receive the archive
    Comm::broadcast(&d_buffer_size, 1, 0);
    Comm::broadcast((char*)d_buffer.data(), d_buffer_size, 0);

    // Setup buffer and archive
    typedef boost::iostreams::basic_array_source<char> source_t;
    source_t device(d_buffer.data(), d_buffer_size);
    boost::iostreams::stream<source_t> s(device);
    boost::archive::binary_iarchive input_archive(s);

    // Fill the node
    input_archive >> nodes;
  }

  // Postconditions
  Ensure(nodes);
}

} // end namespace erme_geometry

//---------------------------------------------------------------------------//
//              end of file NodePartitioner.cc
//---------------------------------------------------------------------------//
