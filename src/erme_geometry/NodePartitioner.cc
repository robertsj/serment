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

  // \todo Make sure we set the proper communicator!

  // Partition the nodes using Comm's default partitioner
  size_t number_nodes = 0;
  size_t local_number_nodes = 0;
  size_t lower_bound = 0;
  if (Comm::rank() == 0) number_nodes = nodes->number_global_nodes();
  Comm::partition(number_nodes, lower_bound, local_number_nodes);

  // Broadcast the nodes.
  broadcast_nodes(nodes);
  Assert(nodes);

  // Set node bounds and finalize.
  nodes->set_bounds(lower_bound, lower_bound + local_number_nodes);
  nodes->finalize();

}

void NodePartitioner::broadcast_nodes(SP_nodelist &nodes)
{

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

}

} // end namespace erme_geometry

//---------------------------------------------------------------------------//
//              end of file NodePartitioner.cc
//---------------------------------------------------------------------------//
