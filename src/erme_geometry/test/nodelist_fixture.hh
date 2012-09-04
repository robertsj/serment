//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   nodelist_fixture.hh
 * \brief  nodelist_fixture 
 * \author Jeremy Roberts
 * \date   Aug 26, 2012
 *
 * Just as each node type should have a prebuilt test fixture, a node
 * list of such types (perhaps heterogeneous) shouls also be created.
 *
 */
//---------------------------------------------------------------------------//

#ifndef NODELIST_FIXTURE_HH_
#define NODELIST_FIXTURE_HH_

#include "erme_geometry/NodeList.hh"
#include "erme_geometry/DummyNode.hh"
#include "erme_geometry/test/node_fixture.hh"

namespace erme_geometry
{

//---------------------------------------------------------------------------//
// DETRAN NODE LISTS
//---------------------------------------------------------------------------//

// Get a list of Detran nodes
NodeList::SP_nodelist cartesian_node_detran_list_2d()
{

  // Create node list
  NodeList::SP_nodelist nodes(new NodeList());

  // Get four, two-dimensional Cartesian test nodes
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
  // | 0 | 1 |   with vacuum on all global surfaces
  //  --- ---
  neigh0[CartesianNode::NORTH]  = NeighborSurface(2, CartesianNode::SOUTH);
  neigh0[CartesianNode::EAST]   = NeighborSurface(1, CartesianNode::WEST);
  neigh1[CartesianNode::WEST]   = NeighborSurface(0, CartesianNode::EAST);
  neigh1[CartesianNode::NORTH]  = NeighborSurface(3, CartesianNode::SOUTH);
  neigh2[CartesianNode::SOUTH]  = NeighborSurface(0, CartesianNode::NORTH);
  neigh2[CartesianNode::EAST]   = NeighborSurface(3, CartesianNode::WEST);
  neigh3[CartesianNode::SOUTH]  = NeighborSurface(1, CartesianNode::NORTH);
  neigh3[CartesianNode::WEST]   = NeighborSurface(2, CartesianNode::EAST);

  // Add nodes
  nodes->add_node(node0, neigh0);
  nodes->add_node(node1, neigh1);
  nodes->add_node(node2, neigh2);
  nodes->add_node(node3, neigh3);
  nodes->finalize();

  return nodes;
}

//---------------------------------------------------------------------------//
// DUMMY NODE LISTS
//---------------------------------------------------------------------------//

// Get a list of Dummy nodes
NodeList::SP_nodelist cartesian_node_dummy_list_1d()
{

  // Create node list
  NodeList::SP_nodelist nodes(new NodeList());;

  // Get four, two-dimensional Cartesian test nodes
  NodeList::SP_node node0(new CartesianNodeDummy(1, 0, 0, 2, 0, 0));
  NodeList::SP_node node1(new CartesianNodeDummy(1, 1, 0, 2, 0, 0));
  NodeList::SP_node node2(new CartesianNodeDummy(1, 2, 0, 2, 0, 0));
  NodeList::SP_node node3(new CartesianNodeDummy(1, 3, 0, 2, 0, 0));

  // Create neighbor lists.
  NodeList::vec_neighbor neigh0(4, NeighborSurface(Node::VACUUM, 0));
  NodeList::vec_neighbor neigh1(4, NeighborSurface(Node::VACUUM, 0));
  NodeList::vec_neighbor neigh2(4, NeighborSurface(Node::VACUUM, 0));
  NodeList::vec_neighbor neigh3(4, NeighborSurface(Node::VACUUM, 0));
  //  --- --- --- ---
  // | 0 | 1 | 2 | 3 | with vacuum on all global surfaces
  //  --- --- --- ---
  neigh0[CartesianNode::EAST] = NeighborSurface(1, CartesianNode::WEST);
  neigh1[CartesianNode::WEST] = NeighborSurface(0, CartesianNode::EAST);

  neigh1[CartesianNode::EAST] = NeighborSurface(2, CartesianNode::WEST);
  neigh2[CartesianNode::WEST] = NeighborSurface(1, CartesianNode::EAST);

  neigh2[CartesianNode::EAST] = NeighborSurface(3, CartesianNode::WEST);
  neigh3[CartesianNode::WEST] = NeighborSurface(2, CartesianNode::EAST);

  // Add nodes
  nodes->add_node(node0, neigh0);
  nodes->add_node(node1, neigh1);
  nodes->add_node(node2, neigh2);
  nodes->add_node(node3, neigh3);
  nodes->finalize();

  return nodes;
}

// Get a list of Dummy nodes
NodeList::SP_nodelist
cartesian_node_dummy_list_2d(int so = 4, int ao = 2, int po = 2)
{

  // Create node list
  NodeList::SP_nodelist nodes(new NodeList());;

  // Get four, two-dimensional Cartesian test nodes
  NodeList::SP_node node0(new CartesianNodeDummy(2, 0, so, po, ao, 0));
  NodeList::SP_node node1(new CartesianNodeDummy(2, 1, so, po, ao, 0));
  NodeList::SP_node node2(new CartesianNodeDummy(2, 2, so, po, ao, 0));
  NodeList::SP_node node3(new CartesianNodeDummy(2, 3, so, po, ao, 0));

  // Create neighbor lists.
  NodeList::vec_neighbor neigh0(4, NeighborSurface(Node::VACUUM, 0));
  NodeList::vec_neighbor neigh1(4, NeighborSurface(Node::VACUUM, 0));
  NodeList::vec_neighbor neigh2(4, NeighborSurface(Node::VACUUM, 0));
  NodeList::vec_neighbor neigh3(4, NeighborSurface(Node::VACUUM, 0));
  //  --- ---
  // | 2 | 3 |
  //  --- ---
  // | 0 | 1 |   with vacuum on all global surfaces
  //  --- ---
  neigh0[CartesianNode::NORTH]  = NeighborSurface(2, CartesianNode::SOUTH);
  neigh0[CartesianNode::EAST]   = NeighborSurface(1, CartesianNode::WEST);
  neigh1[CartesianNode::WEST]   = NeighborSurface(0, CartesianNode::EAST);
  neigh1[CartesianNode::NORTH]  = NeighborSurface(3, CartesianNode::SOUTH);
  neigh2[CartesianNode::SOUTH]  = NeighborSurface(0, CartesianNode::NORTH);
  neigh2[CartesianNode::EAST]   = NeighborSurface(3, CartesianNode::WEST);
  neigh3[CartesianNode::SOUTH]  = NeighborSurface(1, CartesianNode::NORTH);
  neigh3[CartesianNode::WEST]   = NeighborSurface(2, CartesianNode::EAST);

  // Add nodes
  nodes->add_node(node0, neigh0);
  nodes->add_node(node1, neigh1);
  nodes->add_node(node2, neigh2);
  nodes->add_node(node3, neigh3);
  nodes->finalize();

  return nodes;
}

// Get a list of Dummy nodes
NodeList::SP_nodelist cartesian_node_dummy_list_3d()
{

  // Create node list
  NodeList::SP_nodelist nodes(new NodeList());;

  // Get four, two-dimensional Cartesian test nodes
  NodeList::SP_node node0(new CartesianNodeDummy(3, 0, 4, 2, 2, 0));
  NodeList::SP_node node1(new CartesianNodeDummy(3, 1, 4, 2, 2, 0));
  NodeList::SP_node node2(new CartesianNodeDummy(3, 2, 4, 2, 2, 0));
  NodeList::SP_node node3(new CartesianNodeDummy(3, 3, 4, 2, 2, 0));

  // Create neighbor lists.
  NodeList::vec_neighbor neigh0(4, NeighborSurface(Node::VACUUM, 0));
  NodeList::vec_neighbor neigh1(4, NeighborSurface(Node::VACUUM, 0));
  NodeList::vec_neighbor neigh2(4, NeighborSurface(Node::VACUUM, 0));
  NodeList::vec_neighbor neigh3(4, NeighborSurface(Node::VACUUM, 0));
  //  --- ---
  // | 2 | 3 |
  //  --- ---
  // | 0 | 1 |   with vacuum on all global surfaces
  //  --- ---
  neigh0[CartesianNode::NORTH]  = NeighborSurface(2, CartesianNode::SOUTH);
  neigh0[CartesianNode::EAST]   = NeighborSurface(1, CartesianNode::WEST);
  neigh1[CartesianNode::WEST]   = NeighborSurface(0, CartesianNode::EAST);
  neigh1[CartesianNode::NORTH]  = NeighborSurface(3, CartesianNode::SOUTH);
  neigh2[CartesianNode::SOUTH]  = NeighborSurface(0, CartesianNode::NORTH);
  neigh2[CartesianNode::EAST]   = NeighborSurface(3, CartesianNode::WEST);
  neigh3[CartesianNode::SOUTH]  = NeighborSurface(1, CartesianNode::NORTH);
  neigh3[CartesianNode::WEST]   = NeighborSurface(2, CartesianNode::EAST);

  // Add nodes
  nodes->add_node(node0, neigh0);
  nodes->add_node(node1, neigh1);
  nodes->add_node(node2, neigh2);
  nodes->add_node(node3, neigh3);
  nodes->finalize();

  return nodes;
}


} // end namespace erme_geometry

#endif // NODELIST_FIXTURE_HH_ 

//---------------------------------------------------------------------------//
//              end of file nodelist_fixture.hh
//---------------------------------------------------------------------------//
