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

using serment_comm::Comm;

//---------------------------------------------------------------------------//
PostProcess::PostProcess(SP_manager manager)
{
  Require(manager);
  d_server = manager->get_server();
  d_nodes = manager->get_nodes();
  d_indexer = manager->get_indexer();
  d_state = manager->get_state();
  if (Comm::is_global())
    d_F = manager->get_responses()->F;
}

//---------------------------------------------------------------------------//
PostProcess::vec_dbl PostProcess::nodal_fission_density(const double norm)
{
  int N = d_nodes->number_global_nodes();
  vec_dbl fd(N, 0.0);
  double fd_sum = 0.0;

  if (Comm::is_global())
  {
    int i = 0;
    for (int n = d_nodes->lower_bound(); n < d_nodes->upper_bound(); ++n)
    {
      int un = d_nodes->unique_global_index_from_global(n);
      for (int m = 0; m < d_indexer->number_node_moments(un); ++m, ++i)
      {
        fd[n] += (*d_F)[i] * (*d_state->moments())[i];
      }
    }
    for (int n = d_nodes->lower_bound(); n < d_nodes->upper_bound(); ++n)
    {
      fd_sum += fd[n];
    }
  }

  // Gets the complete vector and its sum on all global processes
  Comm::global_sum(fd_sum);
  Comm::global_sum(&fd[0], N);
  for (int n = 0; n < N; ++n)
    fd[n] *= norm / fd_sum;

  return fd;
}

//---------------------------------------------------------------------------//
PostProcess::vec_dbl PostProcess::nodal_power(const double norm)
{
  int N = d_nodes->number_global_nodes();
  vec_dbl power(N, 0.0);
  double power_sum = 0.0;

  if (Comm::is_global())
  {
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
    for (int gn = d_nodes->lower_bound(); gn < d_nodes->upper_bound(); ++gn)
    {
      power_sum += power[gn];
    }
  }
  std::cout << "******" << std::endl;
  serment_comm::Comm::global_sum(power_sum);
  std::cout << "******" << std::endl;

  serment_comm::Comm::global_sum(&power[0], N);
  std::cout << "******" << std::endl;

  for (int n = 0; n < N; ++n)
    power[n] *= norm / power_sum;
  return power;
}

//---------------------------------------------------------------------------//
PostProcess::vec2_dbl PostProcess::pin_power(const double norm)
{
  int N = d_nodes->number_global_nodes();
  vec2_dbl power(N);
  double power_sum = 0.0;


  for (int gn = 0; gn < d_nodes->number_global_nodes(); ++gn)
  {
    power[gn].resize(d_nodes->node(gn)->number_pins(), 0.0);
  }

  if (Comm::is_global())
  {
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
          power_sum += power[gn][p];
        }
      }
    }
  }
  for (int gn = 0; gn < d_nodes->number_global_nodes(); ++gn)
  {
    Comm::global_sum(&power[gn][0], N);
  }
  Comm::global_sum(power_sum);

  for (int gn = 0; gn < d_nodes->number_global_nodes(); ++gn)
    for (int p = 0; p < d_nodes->node(gn)->number_pins(); ++p)
      power[gn][p] *= norm / power_sum;

  return power;
}

} // end namespace erme_utils


