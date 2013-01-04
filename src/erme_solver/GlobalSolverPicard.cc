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
  d_innersolver = new linear_algebra::EigenSolver(d_MR);

  // Check if we are using a non-default keff update
  std::string updater = "default";
  if (d_db->check("picard_update"))
    std::string updater = d_db->get<std::string>("picard_update");
  if (updater == "default")
    d_update = new EigenvalueUpdate();
  else if (updater == "steffensen")
    d_update = new SteffensenUpdate();
  else
    THROW("Unknown Picard eigenvalue updater: " + updater);

  Ensure(d_MR);
  Ensure(d_innersolver);
}

} // end namespace erme_solver

//---------------------------------------------------------------------------//
//              end of file GlobalSolverPicard.cc
//---------------------------------------------------------------------------//
