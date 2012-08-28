//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   nodelist_fixture.hh
 * \brief  nodelist_fixture 
 * \author Jeremy Roberts
 * \date   Aug 26, 2012
 */
//---------------------------------------------------------------------------//

#ifndef NODELIST_FIXTURE_HH_
#define NODELIST_FIXTURE_HH_

#include "erme_geometry/NodeList.hh"
#include "erme_geometry/test/node_fixture.hh"

namespace erme_geometry
{

// Get a list of Detran nodes
NodeList cartesian_node_detran_list_2d()
{

  // Create node list
  NodeList nodes;

  // Get 4 test nodes
  NodeList::SP_node node0 = cartesian_node_detran(2);
  NodeList::SP_node node1 = cartesian_node_detran(2);
  NodeList::SP_node node2 = cartesian_node_detran(2);
  NodeList::SP_node node3 = cartesian_node_detran(2);

  // Create neighbor lists.
  NodeList::vec_neighbor neigh0(4, NeighborSurface(Node::VACUUM, 0));
  NodeList::vec_neighbor neigh1(4, NeighborSurface(Node::VACUUM, 0));
  NodeList::vec_neighbor neigh2(4, NeighborSurface(Node::VACUUM, 0));
  NodeList::vec_neighbor neigh3(4, NeighborSurface(Node::VACUUM, 0));
  //  --- ---
  // | 2 | 3 |
  //  --- ---
  // | 0 | 1 |   with vacuum everywhere else
  //  --- ---
  neigh0[CartesianNode::TOP]    = NeighborSurface(2, CartesianNode::BOTTOM);
  neigh0[CartesianNode::RIGHT]  = NeighborSurface(1, CartesianNode::LEFT);
  neigh1[CartesianNode::LEFT]   = NeighborSurface(0, CartesianNode::RIGHT);
  neigh1[CartesianNode::TOP]    = NeighborSurface(3, CartesianNode::BOTTOM);
  neigh2[CartesianNode::BOTTOM] = NeighborSurface(0, CartesianNode::TOP);
  neigh2[CartesianNode::RIGHT]  = NeighborSurface(3, CartesianNode::RIGHT);
  neigh3[CartesianNode::BOTTOM] = NeighborSurface(1, CartesianNode::TOP);
  neigh3[CartesianNode::LEFT]   = NeighborSurface(2, CartesianNode::RIGHT);

  // Add nodes
  nodes.add_node(node0, neigh0);
  nodes.add_node(node1, neigh1);
  nodes.add_node(node2, neigh2);
  nodes.add_node(node3, neigh3);
  nodes.finalize();

  return nodes;
}

} // end namespace erme_geometry

#endif // NODELIST_FIXTURE_HH_ 

//---------------------------------------------------------------------------//
//              end of file nodelist_fixture.hh
//---------------------------------------------------------------------------//
