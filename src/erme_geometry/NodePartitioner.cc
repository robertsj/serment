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

void NodePartitioner::partition(NodeList &nodes)
{

  // \todo Make sure we set the proper communicator!

  int number_nodes = 0;
  std::vector<int> number_per_process(Comm::size(), 0);

  // Master
  if (Comm::rank() == 0)
  {
    // Preconditions
    Insist(nodes.number_global_nodes() > 0,
           "Node list must contain nodes to partition.");

    number_nodes = nodes.number_global_nodes();

    // Initial guess for nodes per process and the remainder.
    int npp = number_nodes / Comm::size();
    int remainder = number_nodes - npp * Comm::size();

    // Assign the number per process, putting the extras on
    // the processes in reverse.  If there are more processes
    // than nodes, put the nodes on the first processes.
    int bound = Comm::size();
    if (!npp)
    {
      bound = number_nodes;
      npp   = 1;
      remainder = 0;
    }

    for (int i = 0; i < bound; i++)
      number_per_process[i] = npp;

    for (int i = 1; i <= remainder; i++)
      number_per_process[Comm::size()-i] += 1;
  }

  // Broadcast the number of nodes in the problem
  Comm::broadcast(&number_nodes, 1, 0);

  // Broadcast the number of nodes per process.
  Comm::broadcast(&number_per_process[0], number_per_process.size(), 0);

  // Broadcast the nodes
  broadcast_nodes(nodes);

  // Set the appropriate bounds and finalize
  int lb = 0;
  for (int i = 0; i < Comm::rank(); i++)
    lb += number_per_process[i];
  nodes.set_bounds(lb, lb + number_per_process[Comm::rank()]);
  nodes.finalize();

}

void NodePartitioner::broadcast_nodes(NodeList &nodes)
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
