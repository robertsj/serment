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
#define TEST_LIST           \
        FUNC(test_Connect)

#include "TestDriver.hh"
#include "Connect.hh"
#include "ResponseMatrix.hh"
#include "erme_geometry/NodePartitioner.hh"
#include "erme_geometry/test/nodelist_fixture.hh"
#include <iostream>

// Setup

using namespace erme;
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
  typedef serment_comm::Comm Comm;

  Comm::initialize(argc, argv);
  PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

  // Get the node list
  erme_geometry::NodeList::SP_nodelist nodes;
  if (Comm::rank() == 0)
    nodes = erme_geometry::cartesian_node_dummy_list_2d();

  // Partition the nodes
  erme_geometry::NodePartitioner partitioner;
  partitioner.partition(nodes);

  // Create parameter database
  erme_response::ResponseIndexer::SP_db db(new detran::InputDB());
  db->put<int>("dimension", 2);
  db->put<int>("erme_order_reduction", 3);

  // Create indexer
  erme_response::ResponseIndexer::SP_indexer
    indexer(new erme_response::ResponseIndexer(db, nodes));

  // Create new scope for M, since M must be destructed before finalization.
  {
    // Connect
    Connect M(nodes, indexer);

    // Write to binary for inspection in MATLAB
    M.display(Connect::BINARY, "binary.out");
  }

  PetscFinalize();
  Comm::finalize();
  return 0;

}

//---------------------------------------------------------------------------//
//              end of test_StateERME.cc
//---------------------------------------------------------------------------//
