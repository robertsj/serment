//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_MR_scaling.cc
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 * \brief  Tests the parallel scaling of the action M * R
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST             \
        FUNC(test_MR_scaling)

#include "utilities/TestDriver.hh"
#include "Connect.hh"
#include "ResponseMatrix.hh"
#include "erme_geometry/NodePartitioner.hh"
#include "erme_geometry/test/nodelist_fixture.hh"
#include "linear_algebra/LinearAlgebraSetup.hh"
#include <iostream>
#include <boost/program_options.hpp>
#include <cmath>
#include <cstdio>

// Setup

using namespace erme;
using namespace detran_test;
using detran_utilities::soft_equiv;
using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

/*
 *  The global problem has at its core the action M*R, where M
 *  is the connectivity matrix and R is the response matrix.  While
 *  we anticipate the local response computation will represent the
 *  overwhelming majority of all computation, understanding the
 *  scaling of the global problem is important for running the
 *  problem on large machines and for problems in which the
 *  local responses are relatively cheap (e.g. few group diffusion).
 *
 */
int test_MR_scaling(int argc, char *argv[])
{

  //-------------------------------------------------------------------------//
  // INITIALIZE COMM
  //-------------------------------------------------------------------------//

  typedef serment_comm::Comm Comm;
  Comm::initialize(argc, argv);

  //-------------------------------------------------------------------------//
  // GET COMMAND LINE OPTIONS
  //-------------------------------------------------------------------------//

  namespace po = boost::program_options;

  int so;   // spatial order (4 * so = block size)
  int nn;   // number of nodes in one dimension (nn * nn total)
  int nl;   // number of local groups
  int id;   // test id
  int nt;   // number steps for timing loop

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "show this display")
    ("so", po::value< int >(&so)->default_value(0), "spatial order (block size = 4 * so)")
    ("nn", po::value< int >(&nn)->default_value(2), "problem dimension, nn * nn nodes")
    ("nl", po::value< int >(&nl)->default_value(Comm::size()), "number of local groups (must be <= number processes")
    ("id", po::value< int >(&id)->default_value(0), "test id")
    ("nt", po::value< int >(&nt)->default_value(1), "number of timing steps")
  ;

  po::variables_map vm;
  //po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  po::store(po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm);

  po::notify(vm);

  if (vm.count("help"))
  {
    if (Comm::rank() == 0) cout << "Usage: options_description [options]\n";
    if (Comm::rank() == 0) cout << desc;
    Comm::finalize();
    return 0;
  }

  //-------------------------------------------------------------------------//
  // SETUP COMMUNICATORS AND LINEAR ALGEBRA
  //-------------------------------------------------------------------------//

  // Setup local and global communicators.
  Insist(nl <= Comm::size(),
    "Number of local groups must be less than world size");
  Comm::setup_communicators(nl);

  // This sets the PETSc communicator to global
  linear_algebra::initialize(argc, argv);

  //-------------------------------------------------------------------------//
  // SETUP NODES AND PARTITION
  //-------------------------------------------------------------------------//

  // Get the node list
  erme_geometry::NodeList::SP_nodelist nodes;
  if (Comm::rank() == 0)
  {
    //nodes = erme_geometry::cartesian_node_dummy_list_2d(so, 0, 0);
    nodes = erme_geometry::cartesian_node_dummy_list_3d_variable(so, nn);
  }

  // Partition the nodes
  erme_geometry::NodePartitioner partitioner;
  partitioner.partition(nodes);

  TEST(nodes->number_global_nodes() == nn * nn * nn);
  TEST(nodes->number_local_nodes()  >= nn * nn * nn / Comm::size());

  //-------------------------------------------------------------------------//
  // SETUP INDEXER AND SERVER
  //-------------------------------------------------------------------------//

  // Create parameter database
  erme_response::ResponseIndexer::SP_db db(new detran_utilities::InputDB());
  db->put<int>("dimension", 3);
  db->put<int>("erme_order_reduction", 0);

  // Create indexer
  erme_response::ResponseIndexer::SP_indexer
    indexer(new erme_response::ResponseIndexer(db, nodes));

  // Block size
  int bs = indexer->number_node_moments(0);

  // Global problem size
  int gs = indexer->number_global_moments();

  // Create server
  erme_response::ResponseServer::SP_server
    server(new erme_response::ResponseServer(nodes, indexer));

  /*
   *  Update the server, which leads to the first local computations.  In
   *  some cases, this might be the only true computation.  For instance,
   *  if we define the responses as an expansion (R = R0 + R1/k + ...), then
   *  this would compute the expansion coefficients once, store them
   *  locally, and interpolate/expand on-the-fly.
   */
  server->update(1.0);

  //-------------------------------------------------------------------------//
  // BUILD M AND R AND DO TIMING TESTS
  //-------------------------------------------------------------------------//

  std::vector<double> times(nt, 0.0);
  double average_time = 0.0;
  double standard_deviation = 0.0;
  int np = Comm::size();
  // Create R only on the global communicator
  if (Comm::is_global())
  {
    Comm::set(serment_comm::global);

    // Connectivity matrix
    Connect M(nodes, indexer);

    // Response matrix
    ResponseMatrix R(nodes, indexer, server);

    // Update.  Since this uses the dummy data, this simply fills R.
    R.update();

    // Write to binary for inspection in MATLAB
    M.display(ResponseMatrix::BINARY, "connectivity2.out");
    R.display(ResponseMatrix::BINARY, "response_matrix2.out");

    // Create vectors.
    linear_algebra::Vector x(indexer->number_local_moments(), 1.0);
    linear_algebra::Vector y(indexer->number_local_moments(), 0.0);
    linear_algebra::Vector z(indexer->number_local_moments(), 0.0);

    // Do time trials
    for (int trial = 0; trial < nt; trial++)
    {
      Comm::tic();

      // R * x
      R.multiply(x, y);
      // M * y
      M.multiply(y, z);

      times[trial] = Comm::toc();
      average_time += times[trial];
    }
    average_time /= nt;
    for (int trial = 0; trial < nt; trial++)
    {
      standard_deviation +=
        (times[trial] - average_time) * (times[trial] - average_time);
    }
    standard_deviation = std::sqrt(standard_deviation);

    if (Comm::rank() == 0)
    {
      printf("%5i %5i %5i %5i %8i %12.8e  %12.8e \n", np, nl,
        nn*nn*nn, bs, gs, average_time, standard_deviation);
    }

    Comm::set(serment_comm::world);
  }

  //-------------------------------------------------------------------------//
  // WRAP UP
  //-------------------------------------------------------------------------//

  linear_algebra::finalize();
  Comm::finalize();
  return 0;

}

//---------------------------------------------------------------------------//
//              end of test_StateERME.cc
//---------------------------------------------------------------------------//
