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

namespace erme_geometry
{


/*
 *  This implementation views the node list as a one dimensional array
 *  to be distributed.  Hence, each process is given a contiguous
 *  chunk of the global node list.  This is a very easy approach,
 *  but it might not be the best.  Other options would be to
 *  visualize the list as a two or three dimensional array, akin to
 *  the physical model the nodes represent.  Then, processes could
 *  be assigned 2D or 3D contiguous chunks.  That might be a more
 *  efficient approach with respect to memory.
 *
 *  In any case, the partitioning itself should not be expensive,
 *  since the nodes have little data (unless large cross section
 *  sets are hidden within).
 *
 *  \todo Before the partitioning process, the node list should
 *        really be verified.  What this implies is that each
 *        neighbor y of a node x should have node x also listed
 *        as a neighbor.  This is far easier to accomplish by
 *        the master process.  The orders can also be verified,
 *        meaning that the surface between x and y must have
 *        identical moment orders.  If they don't, then some
 *        either error out or adapt.
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
    Insist(nodes.is_finalized(),
           "Node list must be finalized before partitioning");
    Insist(nodes.number_nodes() > 0,
           "Node list must contain nodes to partition.");

    number_nodes = nodes.number_nodes();

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

  // Compute the global offset
  int global_offset = 0;
  for (int i = 0; i < Comm::rank(); i++)
    global_offset += number_per_process[i];

  // Distribute the nodes.
  if (Comm::rank() == 0)
  {
    // Rank 0 keeps its own processes, so start the index
    // at rank 1's nodes
    int global_n = number_per_process[0];

    for (int i = 1; i < Comm::size(); i++)
    {
      for (int n = 0; n < number_per_process[i]; n++, global_n++)
      {
        // Send a node and its neighbors
        send_node(nodes.node(global_n), i);
        Comm::send(&nodes.neighbor(global_n)[0],
                   nodes.neighbor(global_n).size(), i);
      }
    }
  }
  else
  {
    for (int n = 0; n < number_per_process[Comm::rank()]; n++)
    {
      SP_node node;
      // Receive node and its neighbors
      receive_node(node);
      NodeList::vec_int neighbors(node->number_surfaces(), 0);
      Comm::receive(&neighbors[0], 4, 0);
      // Insert
      nodes.add_node(node, neighbors);
    }
  }
  // Rank 0 keeps only its nodes
  if (Comm::rank() == 0)
    nodes.resize(number_per_process[0]);
  nodes.finalize();
}

void NodePartitioner::send_node(SP_node node, int destination)
{
  // Preconditions
  Require(Comm::rank() == 0);
  Require(node);

  // Clear the buffer.  This wipes the contents, but keeps the
  // memory allocated.
  d_buffer.clear();

  // Setup buffer and archive
  typedef boost::iostreams::back_insert_device<std::string> insert_t;
  insert_t inserter(d_buffer);
  boost::iostreams::stream<insert_t> s(inserter);
  boost::archive::binary_oarchive output_archive(s);

  // Archive the node
  output_archive << node;
  s.flush();

  // Send the archive
  int flag = 0;
  d_buffer_size = d_buffer.size();
  flag = Comm::send(&d_buffer_size, 1, destination);
  flag = Comm::send((char*)d_buffer.data(), d_buffer_size, destination);

  // Postconditions
  Require(flag == serment_comm::COMM_SUCCESS);
}

void NodePartitioner::receive_node(SP_node &node)
{
  // Preconditions
  Require(Comm::rank() > 0);
  Require(!node);

  // Clear the buffer
  d_buffer.clear();

  // Receive the archive
  Comm::receive(&d_buffer_size, 1, 0);
  Comm::receive((char*)d_buffer.data(), d_buffer_size, 0);

  // Setup buffer and archive
  typedef boost::iostreams::basic_array_source<char> source_t;
  source_t device(d_buffer.data(), d_buffer_size);
  boost::iostreams::stream<source_t> s(device);
  boost::archive::binary_iarchive input_archive(s);

  // Fill the node
  input_archive >> node;

  // Postconditions
  Ensure(node);
}

} // end namespace erme_geometry

//---------------------------------------------------------------------------//
//              end of file NodePartitioner.cc
//---------------------------------------------------------------------------//
