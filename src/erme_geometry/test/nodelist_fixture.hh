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
  NodeList::SP_nodelist nodes(new NodeList());

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
  NodeList::SP_nodelist nodes(new NodeList());

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
  NodeList::SP_nodelist nodes(new NodeList());

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

// Get a list of Dummy nodes.  This node list provides
// a square map of nodes with varying spatial order.
NodeList::SP_nodelist
cartesian_node_dummy_list_2d_variable(int so, int N)
{

  // Create node list
  NodeList::SP_nodelist nodes(new NodeList());

  // Define one node to use everywhere.
  NodeList::SP_node node(new CartesianNodeDummy(2, 0, so, 0, 0, 0));

  // Vector of neighbor lists.
  std::vector<NodeList::vec_neighbor>
   neighbors(N * N,
             NodeList::vec_neighbor(4, NeighborSurface(Node::VACUUM, 0)));

  //  --- ---
  // | N |N+1| . . .
  //  --- ---
  // | 0 | 1 |...|N-1|
  //  --- --- --- ---

  for (int j = 0; j < N; j++)
  {
    for (int i = 0; i < N; i++)
    {
      int n = i + j * N;
      if (i > 0)   neighbors[n][CartesianNode::WEST] = NeighborSurface(n-1, CartesianNode::EAST);
      if (i < N-1) neighbors[n][CartesianNode::EAST] = NeighborSurface(n+1, CartesianNode::WEST);
      if (j > 0)   neighbors[n][CartesianNode::SOUTH] = NeighborSurface(i + (j-1)*N, CartesianNode::NORTH);
      if (j < N-1) neighbors[n][CartesianNode::NORTH] = NeighborSurface(i + (j+1)*N, CartesianNode::SOUTH);

      nodes->add_node(node, neighbors[n]);
    }
  }
  nodes->finalize();
  return nodes;
}

// Get a list of Dummy nodes.  This node list provides
// a square map of nodes with varying spatial order.
NodeList::SP_nodelist
cartesian_node_dummy_list_3d_variable(int so, int N)
{
  typedef CartesianNode C_N;


  // Create node list
  NodeList::SP_nodelist nodes(new NodeList());

  // Define one node to use everywhere.
  NodeList::SP_node node(new CartesianNodeDummy(3, 0, so, 0, 0, 0));

  // Vector of neighbor lists.
  std::vector<NodeList::vec_neighbor>
   neighbors(N * N * N,
             NodeList::vec_neighbor(6, NeighborSurface(Node::VACUUM, 0)));

  //  --- ---
  // | N |N+1| . . .
  //  --- ---
  // | 0 | 1 |...|N-1|
  //  --- --- --- ---

  for (int k = 0; k < N; k++)
  {
    for (int j = 0; j < N; j++)
    {
      for (int i = 0; i < N; i++)
      {
        // cardinal index
        int me = i + j * N + k * N * N;
        int w = me - 1;
        int e = me + 1;
        int s = i + (j-1) * N + k * N * N;
        int n = i + (j+1) * N + k * N * N;
        int b = i + j * N + (k-1) * N * N;
        int t = i + j * N + (k+1) * N * N;

        if (i > 0)     neighbors[me][C_N::WEST]    = NeighborSurface(w, C_N::EAST);
        if (i < N - 1) neighbors[me][C_N::EAST]    = NeighborSurface(e, C_N::WEST);
        if (j > 0)     neighbors[me][C_N::SOUTH]   = NeighborSurface(s, C_N::NORTH);
        if (j < N - 1) neighbors[me][C_N::NORTH]   = NeighborSurface(n, C_N::SOUTH);
        if (k > 0)     neighbors[me][C_N::BOTTOM]  = NeighborSurface(b, C_N::NORTH);
        if (k < N - 1) neighbors[me][C_N::TOP]     = NeighborSurface(t, C_N::BOTTOM);
        nodes->add_node(node, neighbors[me]);
      }
    }
  }
  nodes->finalize();
  return nodes;
}


} // end namespace erme_geometry

#endif // NODELIST_FIXTURE_HH_ 

//---------------------------------------------------------------------------//
//              end of file nodelist_fixture.hh
//---------------------------------------------------------------------------//
