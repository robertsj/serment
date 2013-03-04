//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   GlobalSolverPicard.cc
 *  @brief  GlobalSolverPicard
 *  @author Jeremy Roberts
 *  @date   Oct 1, 2012
 */
//---------------------------------------------------------------------------//

#include "GlobalSolverPicard.hh"
#include "SteffensenUpdate.hh"

namespace erme_solver
{

//---------------------------------------------------------------------------//
GlobalSolverPicard::GlobalSolverPicard(SP_db      db,
                                       SP_indexer indexer,
                                       SP_server  server,
                                       SP_state   state,
                                       SP_R       R,
                                       SP_M       M,
                                       SP_F       F,
                                       SP_A       A,
                                       SP_L       L)
  : Base(db, indexer, server, state, R, M, F, A, L)
  , d_J0(new linear_algebra::Vector(d_state->local_size(), 0.0))
  , d_J1(new linear_algebra::Vector(d_state->local_size(), 0.0))
{
  // Create MR operator
  d_MR = new OperatorMR(d_R, d_M);

  // Create inner eigensolver
  int inner_max_iters = 1e5;
  if (d_db->check("erme_inner_max_iters"))
    inner_max_iters = d_db->get<int>("erme_inner_max_iters");
  double inner_tolerance = 1e-10;
  if (d_db->check("erme_inner_tolerance"))
    inner_tolerance = d_db->get<double>("erme_inner_tolerance");
  d_innersolver = new linear_algebra::EigenSolver(d_MR,
                                                  linear_algebra::Matrix::SP_matrix(0),
                                                  inner_max_iters, inner_tolerance);

  // Check if we are using a non-default keff update
  std::string updater = "default";
  if (d_db->check("picard_update"))
    updater = d_db->get<std::string>("picard_update");
  if (updater == "default")
    d_update = new EigenvalueUpdate();
  else if (updater == "steffensen")
    d_update = new SteffensenUpdate();
  else
    THROW("Unknown Picard eigenvalue updater: " + updater);
  std::cout << "Using eigenvalue updater: " << updater << std::endl;
  // Post conditions
  Ensure(d_MR);
  Ensure(d_innersolver);
}

} // end namespace erme_solver

//---------------------------------------------------------------------------//
//              end of file GlobalSolverPicard.cc
//---------------------------------------------------------------------------//
