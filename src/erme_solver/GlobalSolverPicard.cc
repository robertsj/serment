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
#include <cstdio>

namespace erme_solver
{

//---------------------------------------------------------------------------//
GlobalSolverPicard::GlobalSolverPicard(SP_db      					db,
                                       SP_indexer 					indexer,
                                       SP_server  					server,
                                       SP_state   					state,
                                       SP_responsecontainer responses)
  : Base(db, indexer, server, state, responses)
{

	if (serment_comm::Comm::is_global())
	{
		using namespace linear_algebra;

	  d_J0 = new Vector(d_state->local_size(), 0.0);
	  d_J1 = new Vector(d_state->local_size(), 0.0);

		// Create MR operator
		d_MR = new OperatorMR(d_R, d_M);

		// Inner eigensolver iterations
		int inner_max_iters = 1e5;
		if (d_db->check("erme_inner_max_iters"))
			inner_max_iters = d_db->get<int>("erme_inner_max_iters");

		// Inner eigensolver tolderance
		double inner_tolerance = 1e-10;
		if (d_db->check("erme_inner_tolerance"))
			inner_tolerance = d_db->get<double>("erme_inner_tolerance");

		d_innersolver = new EigenSolver(d_MR,
				                            Matrix::SP_matrix(0),
				                            inner_max_iters,
				                            inner_tolerance);

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
		if (serment_comm::Comm::world_rank() == 0)
		  std::cout << "Using eigenvalue updater: " << updater << std::endl;
	}

}

//---------------------------------------------------------------------------//
void GlobalSolverPicard::solve()
{
	using serment_comm::Comm;
	using std::cout;
	using std::endl;
  Insist(serment_comm::communicator == serment_comm::world,
  		   "Solvers must be accessed on the world communicator.");

  if (Comm::world_rank() == 0) cout << " Picard: solving... " << endl;

  // Get initial keff guess
  double keff = 1.0;
  if (d_db->check("erme_initial_keff"))
    keff = d_db->get<double>("erme_initial_keff");

  // Start lambda at unity
  double lambda = 1.0;

  if (Comm::is_global())
  {
		// Set space-angle zeroth order moments to uniform guess.  Not
  	// used by SLEPc but will in a hand-coded power iteration.
		for (int i = 0; i < d_indexer->number_local_moments(); ++i)
		{
			erme_response::ResponseIndex ri = d_indexer->response_index_from_local(i);
			if (ri.azimuth + ri.polar + ri.space0 + ri.space1 == 0) (*d_J0)[i] = 1.0;
		}
		// Ensure a normalized initial guess and initialize the responses
		d_J0->scale(1.0/d_J0->norm(d_J0->L2));
  }

  // Compute initial responses
  update_response(keff);

  if (Comm::is_global())
  {
		d_R->display(d_R->MATLAB, "R_i.out");
		d_M->display(d_R->MATLAB, "M.out");
		d_L->display(d_L->MATLAB, "L_i.out");
//		d_F->display();
//		d_A->display();
  }

  // Initialize balance parameters
  double loss       = 0;
  double gain       = 0;
  double absorption = 0;
  double leakage    = 0;

  // Count total iterations
  int innertot = 0;

  // Compute the initial residual norm
  double norm = 1.0;
  if (Comm::is_global())
  	norm = d_residual->compute_norm(*d_J0, keff, lambda);
  Comm::broadcast(&norm, 1, 0);
  d_residual_norms.push_back(norm);

  //-------------------------------------------------------------------------//
  // OUTER ITERATIONS
  //-------------------------------------------------------------------------//

  int it = 1; // count outer iteration
  for (; it <= d_maximum_iterations; it++)
  {

  	if (Comm::is_global())
  	{

			//---------------------------------------------------------------------//
			// INNER ITERATIONS -- solves M*R*X = lambda*X
			//---------------------------------------------------------------------//

			lambda = d_innersolver->solve(d_J0);

			//---------------------------------------------------------------------//
			// EIGENVALUE UPDATE -- k = fission / (absorption + leakage)
			//---------------------------------------------------------------------//

			// Compute gains and losses.
			gain        = d_F->dot(*d_J0);
			absorption  = d_A->dot(*d_J0);
			leakage     = d_L->leakage(*d_J0);
			loss        = absorption + leakage;


		  if (Comm::is_global())
		  {
				d_R->display(d_R->MATLAB, "R.out");
				d_L->display(d_L->MATLAB, "L.out");
		//		d_F->display();
		//		d_A->display();
		  }
//			if (serment_comm::Comm::world_rank() == 0)
//			{
//				std::cout << " GAIN = "     << gain << std::endl;
//				std::cout << " ABS = "      << absorption << std::endl;
//				std::cout << " leakage = "  << leakage << std::endl;
//			}

			// Initial update of keff
			keff = gain / loss;

			// Improved the update, if applicable
			keff = d_update->compute(keff, d_J0);

  	} // global
  	Comm::broadcast(&keff, 1, 0);

    // Update the responses
    update_response(keff);

    // Compute the norm of the nonlinear residual.
    if (Comm::is_global())
    	norm = d_residual->compute_norm(*d_J0, keff, lambda);
    Comm::broadcast(&norm, 1, 0);
    d_residual_norms.push_back(norm);

    if (serment_comm::Comm::world_rank() == 0)
    {
    	printf(" PICARD IT: %3i NORM: %8.6e LAMBDA: %12.9f KEFF: %12.9f INNERS %8i \n",
        it, norm, lambda, keff, innertot);
    }

    // Check convergence
    if (norm < d_tolerance) break;

  } // end outers

  if (Comm::world_rank() == 0)
  {
		std::printf(" FINAL PICARD ITERATION EIGENVALUES: \n");
		std::printf(" **** FINAL KEFF        = %12.9f \n", keff);
		std::printf(" **** FINAL LAMBDA      = %12.9f \n", lambda);
		std::printf(" **** OUTER ITERATIONS  = %8i \n", it);
		std::printf(" **** INNER ITERATIONS  = %8i \n", innertot);
  }

  // Update the state
  d_state->update(d_J0, keff, lambda);

  return;
}

} // end namespace erme_solver

//---------------------------------------------------------------------------//
//              end of file GlobalSolverPicard.cc
//---------------------------------------------------------------------------//
