//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   TwoGroupKernel.cc
 * \author Jeremy Roberts
 * \date   10/14/2011
 * \brief  TwoGroupKernel class member definitions.
 * \note   Copyright (C) 2011 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 140                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-09-14 12:53:40 -0400 (Wed, 14 Sep 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#include "Definitions.hh"
#include "TwoGroupKernel.hh"

#include "utilities/Constants.hh"
#include "utilities/DBC.hh"

namespace erme
{

TwoGroupKernel::TwoGroupKernel()
{
  // nothing here for now.
}

TwoGroupKernel::~TwoGroupKernel()
{
  // nothing here for now.
}


void TwoGroupKernel::response(SP_responsefunction responses, Vec_Dbl &data)
{

  // Compute all feasible states, i.e. those with est. keff < 1.0, fill
  // a table for me, and clear the rest.
  Mat_Dbl data_ranges_full;
  Mat_Dbl data_ranges;
  fill_data(data_ranges_full, data_ranges);

  // Initialize all the problems.  Each process has one.
  d_input   = new Diff2dInput();
  d_problem = new Diff2dProblem( *d_input, 0 );
  d_solver  = new Diff2dSolver( *d_input, *d_problem, 0 );


  // Timing.
  double t_total  = 0;
  double t_begin  = 0;
  double t_state  = 0;
  double t_avg    = 0;

  int batch_size = d_up - d_low;
  for (int state = 0; state < batch_size; state++)
  {

    if (d_rank == 0)
      t_begin = MPI_Wtime();

    Vec_Dbl &data = data_ranges[state];
    update_input(data);
    d_problem->setK();

    if (keff(data) >= 1.0)
    {
      // do something to indicate.
    }

    for (int order = 0; order <= d_order; order++)
    {
      for (int group = 0; group < 2; group++)
      {
        // Solve for the flux and responses.
        d_solver->solve();
       // d_solver->compute_responses();
      }
    }

    // Put input and associated responses in output.


    if (d_rank == 0)
    {
      t_state = MPI_Wtime() - t_begin;
      t_total += t_state;
      if ( state % 100 == 0 )
      {
        t_avg = t_total / state;
        printf ("                  Finished: %8i / %8i     \n ", state, 1);
        printf ("          Total time [sec]: %13.10f       \n ", t_total );
        printf ("  Average time/state [sec]: %13.10f       \n ", t_avg   );
        printf ("  Est. time remaining [hr]: %13.10f       \n ",
                t_avg*(batch_size-state)/3600);
      }
    }

  }


}

} // end namespace erme

//---------------------------------------------------------------------------//
//                 end of TwoGroupKernel.hh
//---------------------------------------------------------------------------//

