//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   PostProcess.cc
 *  @author robertsj
 *  @date   Jan 22, 2013
 *  @brief  PostProcess class definition.
 */
//---------------------------------------------------------------------------//

#include "PostProcess.hh"

namespace erme_utils
{

//---------------------------------------------------------------------------//
PostProcess::PostProcess(SP_db db,
                         SP_nodelist nodes,
                         SP_indexer indexer,
                         SP_server server,
                         SP_state state,
                         SP_R R,
                         SP_M M,
                         SP_F F,
                         SP_A A,
                         SP_L L)
  : d_db(db)
  , d_nodes(nodes)
  , d_indexer(indexer)
  , d_server(server)
  , d_state(state)
  , d_R(R)
  , d_M(M)
  , d_F(F)
  , d_A(A)
  , d_L(L)
{

}

//---------------------------------------------------------------------------//
PostProcess::vec_dbl PostProcess::nodal_fission_density(const double norm)
{
  int N = d_nodes->number_global_nodes();
  vec_dbl fd(N, 0.0);

  int i = 0;
  for (int n = d_nodes->lower_bound(); n < d_nodes->upper_bound(); ++n)
  {
    int un = d_nodes->unique_global_index_from_global(n);
    for (int m = 0; m < d_indexer->number_node_moments(un); ++m, ++i)
    {
      fd[n] += (*d_F)[i] * d_state->moments()[i];
    }
  }
  serment_comm::Comm::global_sum(&fd[0], N);
  double fd_sum = 0.0;
  for (int n = 0; n < fd.size(); ++n)
    fd_sum += fd[n];
  for (int n = 0; n < fd.size(); ++n)
    fd[n] *= norm / fd_sum;
  return fd;
}

} // end namespace erme_utils


