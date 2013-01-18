//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_ResponseSourceDetran.cc
 *  @author Jeremy Roberts
 *  @date   Aug 19, 2012
 *  @brief  Test of ResponseServer class.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                               \
        FUNC(test_ResponseSourceDetran_1D_DIFF) \
        FUNC(test_ResponseSourceDetran_2D_DIFF)

#include "utilities/TestDriver.hh"
#include "erme_response/ResponseSourceDetran.hh"
#include "erme_response/ResponseIndexer.hh"
#include "erme_response/ResponseServer.hh"
#include "erme_geometry/NodePartitioner.hh"
#include "erme_geometry/CartesianNodeDetran.hh"
#include "linear_algebra/LinearAlgebraSetup.hh"
#include "geometry/detran_geometry.hh"
#include "boundary/BoundaryDiffusion.hh"
#include "boundary/BoundarySN.hh"
#include <iostream>

// Setup
#include "erme_geometry/test/nodelist_fixture.hh"

using namespace erme_geometry;
using namespace erme_response;
using namespace detran_test;
using namespace detran_utilities;
using namespace detran;
using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//---------------------------------------------------------------------------//
// TEST DEFINITIONS
//---------------------------------------------------------------------------//

int test_ResponseSourceDetran_1D_DIFF(int argc, char *argv[])
{

  //-------------------------------------------------------------------------//
  // SETUP COMM
  //-------------------------------------------------------------------------//

  typedef serment_comm::Comm Comm;
  typedef ResponseSourceDetran<detran::_1D> Source_T;

  Comm::initialize(argc, argv);
  if (Comm::size() > 1) return 0;
  int number_local_comm = 1;
  Comm::setup_communicators(number_local_comm);

  linear_algebra::initialize(argc, argv);

  //-------------------------------------------------------------------------//
  // SETUP DETRAN NODE
  //-------------------------------------------------------------------------//

  // Input
  InputDB::SP_input db(new InputDB());
  db->put<int>("number_groups", 2);
  db->put<int>("dimension", 1);
  db->put<std::string>("equation", "diffusion");

  // Mesh
  vec_int fm(1, 10);
  vec_dbl cm(2, 0);
  cm[1] = 10.0;
  vec_int mt(1, 0);
  detran_geometry::Mesh::SP_mesh mesh;
  mesh = detran_geometry::Mesh1D::Create(fm, cm, mt);

  // Material
  detran_material::Material::SP_material mat;
  mat = detran_material::Material::Create(1, 2);
  mat->set_sigma_t(0, 0, 1.0);
  mat->set_sigma_t(0, 1, 1.0);
  mat->set_sigma_s(0, 1, 0, 0.5);
  mat->set_sigma_s(0, 0, 1, 0.5);
  mat->set_sigma_f(0, 0, 0.5);
  mat->set_chi(0, 0,     1.0);
  mat->compute_diff_coef();
  mat->compute_sigma_a();
  mat->finalize();

  // Node
  vec2_size_t so(2, vec_size_t(0));
  vec_size_t  ao(2, 0);
  vec_size_t  po(2, 0);
  vec_size_t  eo(2, 1);
  vec_dbl     width(3, 1.0); width[0] = 10.0;
  typename Source_T::SP_node node(new erme_geometry::
    CartesianNodeDetran(1, "lala", so, ao, po, eo, width, db, mat, mesh));

  //-------------------------------------------------------------------------//
  // NODELIST
  //-------------------------------------------------------------------------//

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

  //-------------------------------------------------------------------------//
  // RESPONSE AND SOURCE
  //-------------------------------------------------------------------------//

  // Indexer
  db->put<int>("erme_order_reduction", 0);
  ResponseIndexer::SP_indexer indexer(new ResponseIndexer(db, nodes));

  indexer->display();

  // Dummy response -- 2 surfaces x 2 groups
  detran_utilities::size_t N = 2 * 2;
  NodeResponse::SP_response response(new NodeResponse(N, 2));

  // Source
  Source_T source(node, indexer);

  // Solver
  typename Source_T::SP_solver solver = source.solver();
  BoundaryDiffusion<_1D>::SP_boundary b = solver->boundary();
  BoundaryDiffusion<_1D> &boundary = *b;

  source.update(1.0);

  for (int i = 0; i < indexer->number_node_moments(0); ++i)
  {
    ResponseIndex index = indexer->response_index_from_unique_local(i);
    std::cout << " INDEX = " << index << std::endl;
    source.compute(response, index);
  }
  response->display();
  return 0;
}

int test_ResponseSourceDetran_2D_DIFF(int argc, char *argv[])
{

  //-------------------------------------------------------------------------//
  // SETUP COMM
  //-------------------------------------------------------------------------//

  typedef serment_comm::Comm Comm;
  typedef ResponseSourceDetran<detran::_2D> Source_T;

  Comm::initialize(argc, argv);
  if (Comm::size() > 1) return 0;
  int number_local_comm = 1;
  Comm::setup_communicators(number_local_comm);

  linear_algebra::initialize(argc, argv);

  //-------------------------------------------------------------------------//
  // SETUP DETRAN NODE
  //-------------------------------------------------------------------------//

  // Input
  InputDB::SP_input db(new InputDB());
  db->put<int>("number_groups", 2);
  db->put<int>("dimension", 2);
  db->put<std::string>("equation", "diffusion");

  // Mesh
  vec_int fm(1, 10);
  vec_dbl cm(2, 0);
  cm[1] = 10.0;
  vec_int mt(1, 0);
  detran_geometry::Mesh::SP_mesh mesh;
  mesh = detran_geometry::Mesh2D::Create(fm, fm, cm, cm, mt);
  mesh->display();
  // Material
  detran_material::Material::SP_material mat;
  mat = detran_material::Material::Create(1, 2);
  mat->set_sigma_t(0, 0, 1.0);
  mat->set_sigma_t(0, 1, 1.0);
  mat->set_sigma_s(0, 1, 0, 0.5);
  mat->set_sigma_s(0, 0, 1, 0.5);
  mat->set_sigma_f(0, 0, 0.5);
  mat->set_chi(0, 0,     1.0);
  mat->compute_diff_coef();
  mat->compute_sigma_a();
  mat->finalize();

  // Node
  vec2_size_t so(4, vec_size_t(1, 1));
  vec_size_t  ao(4, 0);
  vec_size_t  po(4, 0);
  vec_size_t  eo(4, 1);
  vec_dbl     width(3, 10.0); width[2] = 1.0;
  typename Source_T::SP_node node(new erme_geometry::
    CartesianNodeDetran(2, "lala", so, ao, po, eo, width, db, mat, mesh));

  //-------------------------------------------------------------------------//
  // NODELIST
  //-------------------------------------------------------------------------//

  NodeList::vec2_neighbor
   neighbors(3, NodeList::vec_neighbor(4, NeighborSurface(Node::VACUUM, 0)));
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

  //-------------------------------------------------------------------------//
  // RESPONSE AND SOURCE
  //-------------------------------------------------------------------------//

  // Indexer
  db->put<int>("erme_order_reduction", 0);
  ResponseIndexer::SP_indexer indexer(new ResponseIndexer(db, nodes));

  indexer->display();

  // Dummy response -- 4 surfaces x 2 spatial moments x 2 groups
  detran_utilities::size_t N = 4 * 2 * 2;
  NodeResponse::SP_response response(new NodeResponse(N, 4));

  // Source
  Source_T source(node, indexer);

  // Solver
  typename Source_T::SP_solver solver = source.solver();
  BoundaryDiffusion<_2D>::SP_boundary b = solver->boundary();
  BoundaryDiffusion<_2D> &boundary = *b;

  source.update(1.0);

  for (int i = 0; i < indexer->number_node_moments(0); ++i)
  {
    ResponseIndex index = indexer->response_index_from_unique_local(i);
    std::cout << " INDEX = " << index << std::endl;
    source.compute(response, index);
  }
  response->display();
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_ResponseSourceDetran.cc
//---------------------------------------------------------------------------//
