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
PostProcess::PostProcess(SP_manager manager)
{
  Require(manager);
//  : d_db(db)
//  , d_nodes(nodes)
//  , d_indexer(indexer)
//  , d_server(server)
//  , d_state(state)
//  , d_R(R)
//  , d_M(M)
//  , d_F(F)
//  , d_A(A)
//  , d_L(L)
  d_server = manager->get_server();
  d_nodes = manager->get_nodes();
  d_indexer = manager->get_indexer();
  d_state = manager->get_state();
  d_F = manager->get_responses()->F;
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
      fd[n] += (*d_F)[i] * (*d_state->moments())[i];
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

//---------------------------------------------------------------------------//
PostProcess::vec_dbl PostProcess::nodal_power(const double norm)
{
  int N = d_nodes->number_global_nodes();
  vec_dbl power(N, 0.0);

  int i = 0;
  for (int gn = d_nodes->lower_bound(); gn < d_nodes->upper_bound(); ++gn)
  {
    int ugn = d_nodes->unique_global_index_from_global(gn);
    int  ln = d_nodes->local_index_from_global(gn);
    erme_response::ResponseServer::SP_response rf = d_server->response(ln);
    for (int m = 0; m < d_indexer->number_node_moments(ugn); ++m, ++i)
    {
      power[gn] += rf->nodal_power(m) * (*d_state->moments())[i];
    }
  }
  serment_comm::Comm::global_sum(&power[0], N);
  double power_sum = 0.0;
  for (int n = 0; n < power.size(); ++n)
  {
    power_sum += power[n];
  }
  for (int n = 0; n < power.size(); ++n)
    power[n] *= norm / power_sum;
  return power;
}

//---------------------------------------------------------------------------//
PostProcess::vec2_dbl PostProcess::pin_power(const double norm)
{
  int N = d_nodes->number_global_nodes();
  vec2_dbl power(N);
  for (int gn = d_nodes->lower_bound(); gn < d_nodes->upper_bound(); ++gn)
  {
    power[gn].resize(d_nodes->node(gn)->number_pins(), 0.0);
  }

  int i = 0;
  for (int gn = d_nodes->lower_bound(); gn < d_nodes->upper_bound(); ++gn)
  {
    int ugn = d_nodes->unique_global_index_from_global(gn);
    int  ln = d_nodes->local_index_from_global(gn);
    erme_response::ResponseServer::SP_response rf = d_server->response(ln);
    for (int m = 0; m < d_indexer->number_node_moments(ugn); ++m, ++i)
    {
      for (int p = 0; p < d_nodes->node(gn)->number_pins(); ++p)
      {
        power[gn][p] += rf->pin_power(p, m) * (*d_state->moments())[i];
      }
    }
  }

  //serment_comm::Comm::global_sum(&power[0], N);
  double power_sum = 0.0;
  for (int gn = d_nodes->lower_bound(); gn < d_nodes->upper_bound(); ++gn)
    for (int p = 0; p < d_nodes->node(gn)->number_pins(); ++p)
      power_sum += power[gn][p];
  for (int gn = d_nodes->lower_bound(); gn < d_nodes->upper_bound(); ++gn)
    for (int p = 0; p < d_nodes->node(gn)->number_pins(); ++p)
      power[gn][p] *= norm / power_sum;

  return power;
}

} // end namespace erme_utils


