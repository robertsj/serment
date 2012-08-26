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
  NodeList::vec_int neigh0(4, Node::VACUUM);
  NodeList::vec_int neigh1(4, Node::VACUUM);
  NodeList::vec_int neigh2(4, Node::VACUUM);
  NodeList::vec_int neigh3(4, Node::VACUUM);
  //  --- ---
  // | 2 | 3 |
  //  --- ---
  // | 0 | 1 |   with vacuum everywhere else
  //  --- ---
  neigh0[CartesianNode::TOP]    = 2;
  neigh0[CartesianNode::RIGHT]  = 1;
  neigh1[CartesianNode::LEFT]   = 0;
  neigh1[CartesianNode::TOP]    = 3;
  neigh2[CartesianNode::BOTTOM] = 0;
  neigh2[CartesianNode::RIGHT]  = 3;
  neigh3[CartesianNode::BOTTOM] = 1;
  neigh3[CartesianNode::LEFT]   = 2;

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
