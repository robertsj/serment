//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_NodeList.cc
 *  @author Jeremy Roberts
 *  @date   Aug 19, 2012
 *  @brief  Test of NodeList class.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST         \
        FUNC(test_NodeList)

#include "utilities/TestDriver.hh"
#include "erme_geometry/NodePartitioner.hh"
#include "erme_geometry/NodeList.hh"
#include "erme_geometry/CartesianNode.hh"
#include "erme_geometry/test/nodelist_fixture.hh"
#include <iostream>

using namespace erme_geometry;
using namespace detran_test;
using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

// Test of basic public interface
int test_NodeList(int argc, char *argv[])
{
  // Single process only
  if (!(serment_comm::Comm::rank() == 0)) return 0;

  // Test 1D Cartesian
  {
    // Create 1 Cartesian node.
    Node::vec2_size_t so(2, Node::vec_size_t(1, 0));
    Node::vec_size_t ao(2, 0);
    Node::vec_size_t po(2, 0);
    Node::vec_size_t eo(2, 0);
    Node::vec_dbl w(3, 1.0);
    NodeList::SP_node node(new CartesianNode(1, "testnode", so, ao, po, eo, w));

    // Create neighbor list. Default to vacuum, and then change the
    // 2 shared boundaries.
    NodeList::vec2_neighbor
      neighbors(3, NodeList::vec_neighbor(2, NeighborSurface(Node::VACUUM, 0)));
    neighbors[0][CartesianNode::EAST] = NeighborSurface(1, CartesianNode::WEST);
    neighbors[1][CartesianNode::WEST] = NeighborSurface(0, CartesianNode::EAST);
    neighbors[1][CartesianNode::EAST] = NeighborSurface(2, CartesianNode::WEST);
    neighbors[2][CartesianNode::WEST] = NeighborSurface(1, CartesianNode::EAST);

    // Create a map [0, 0, 0]
    NodeList::vec_int nodemap(3, 0);

    // Node list
    NodeList::SP_nodelist nodes = NodeList::Create();
    nodes->add_node(node);
    nodes->set_nodal_map(nodemap, neighbors);

    // Partition
    NodePartitioner P;
    P.partition(nodes);

    nodes->display();
    cout << " number global nodes = " << nodes->number_global_nodes() << endl;
    cout << " number local nodes  = " << nodes->number_local_nodes() << endl;
    cout << " number global unique nodes = " << nodes->number_unique_global_nodes() << endl;
    cout << " number local unique nodes = " << nodes->number_unique_local_nodes() << endl;

    //
    TEST(nodes->number_global_nodes() == 3);
    TEST(nodes->number_local_nodes()  == 3);
    TEST(nodes->number_unique_global_nodes()  == 1);
    TEST(nodes->number_unique_local_nodes()   == 1);

//    TEST(nodes->number_global_nodes() == 3);
//    TEST(nodes->neighbor(0, CartesianNode::NORTH).neighbor() == 2);
//    TEST(nodes->neighbor(0, CartesianNode::NORTH).surface()  == CartesianNode::SOUTH);
//    TEST(nodes->neighbor(1, CartesianNode::EAST).neighbor()  == Node::VACUUM);
//    TEST(nodes->node(2)->dimension() == 2);
  }
  // Test 2D Cartesian
  {
    NodeList::SP_nodelist nodes = cartesian_node_dummy_list_2d();
    TEST(nodes->number_global_nodes() == 4);
    TEST(nodes->neighbor(0, CartesianNode::NORTH).neighbor() == 2);
    TEST(nodes->neighbor(0, CartesianNode::NORTH).surface()  == CartesianNode::SOUTH);
    TEST(nodes->neighbor(1, CartesianNode::EAST).neighbor()  == Node::VACUUM);
    TEST(nodes->node(2)->dimension() == 2);

    TEST(nodes->lower_bound() == 0);
    TEST(nodes->upper_bound() == 4);

    // there are four unique nodes in the list, but two are indexed
    //   global:  2  2  0  0
    //    local:  2  2  0  0
    //       ug:  0  2
    //       ul:  0  2

    TEST(nodes->local_index_from_global(0) == 0);
    TEST(nodes->local_index_from_global(1) == 1);
    TEST(nodes->local_index_from_global(2) == 2);
    TEST(nodes->local_index_from_global(3) == 3);

    TEST(nodes->global_index_from_local(0) == 0);
    TEST(nodes->global_index_from_local(1) == 1);
    TEST(nodes->global_index_from_local(2) == 2);
    TEST(nodes->global_index_from_local(3) == 3);

    TEST(nodes->unique_global_index_from_global(0) == 2);
    TEST(nodes->unique_global_index_from_global(1) == 2);
    TEST(nodes->unique_global_index_from_global(2) == 0);
    TEST(nodes->unique_global_index_from_global(3) == 0);

    TEST(nodes->unique_global_index_from_unique_local(0) == 0);
    TEST(nodes->unique_global_index_from_unique_local(1) == 2);

    TEST(nodes->unique_local_index_from_unique_global(0) == 0);
    TEST(nodes->unique_local_index_from_unique_global(2) == 1);

    for (int i = 0; i < 4; ++i)
    {
      cout << "  g=" << i
           << "  l=" << nodes->local_index_from_global(i)
           << " ug=" << nodes->unique_global_index_from_global(i)
           << " ul=" << nodes->unique_local_index_from_unique_global(nodes->unique_global_index_from_global(i))
           << endl;
    }


  }
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_NodeList.cc
//---------------------------------------------------------------------------//
