//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ResponseServer.cc
 *  @brief ResponseServer class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "ResponseServer.hh"
#include "ResponseSourceFactory.hh"
#include "comm/Comm.hh"
#include "utilities/SoftEquivalence.hh"
#include <algorithm>
#include <iostream>
#include <cstdio>

namespace erme_response
{

using serment_comm::Comm;
using detran_utilities::soft_equiv;
using std::cout;
using std::endl;


//----------------------------------------------------------------------------//
ResponseServer::ResponseServer(SP_nodelist  nodes,
                               SP_indexer   indexer,
                               std::string  dbname,
                               size_t       dborder,
                               std::string  dbtype)
  : d_nodes(nodes)
  , d_indexer(indexer)
  , d_keff(-1.0)
  , d_keff_1(-2.0)
  , d_is_updated(false)
  , d_number_unique_evaluations(0)
  , d_number_swapped_evaluations(0)
  , d_response_time(0.0)
{
  Require(d_nodes);
  Require(d_indexer);

  d_sources.resize(d_nodes->number_unique_local_nodes());
  d_responses.resize(d_nodes->number_unique_local_nodes());
  d_responses_1.resize(d_nodes->number_unique_local_nodes());

  // Build the response database if requested
  if (dbname != "")
  {
    d_rfdb = ResponseDatabase::Create(dbname, dborder, dbtype);
  }

  /*
   *  Create the sources and responses.  Note, these are
   *  done for *unique* nodes.  Thus, repeated nodes are
   *  computed only once locally.
   */
  ResponseSourceFactory builder;
  for (size_t node_ul = 0; node_ul < d_sources.size(); ++node_ul)
  {
    // Build the sources
    size_t node_ug =  d_nodes->unique_global_index_from_unique_local(node_ul);
    d_sources[node_ul] = builder.build(nodes->unique_node(node_ug), d_indexer);
    Assert(d_sources[node_ul]);

    // Build the nodal response containers
    d_responses[node_ul] =
        new NodeResponse(d_indexer->number_node_moments(node_ug),
                         d_nodes->unique_node(node_ug)->number_surfaces(),
                         d_nodes->unique_node(node_ug)->number_pins());

    d_responses_1[node_ul] =
        new NodeResponse(d_indexer->number_node_moments(node_ug),
                         d_nodes->unique_node(node_ug)->number_surfaces(),
                         d_nodes->unique_node(node_ug)->number_pins());
  }

}

//----------------------------------------------------------------------------//
bool ResponseServer::update(const double keff, int msg)
{
  Require(serment_comm::communicator == serment_comm::world);

  bool store_k = true;

  double et = Comm::wtime();

  Comm::tic();

  // Switch to local communicator.
  Comm::set(serment_comm::local);

  // Broadcast new keff to local communicator.  The local root
  // must have it from the global solve.
  double k = keff;
  Comm::broadcast(&k, 1, 0);

//  std::printf(" UPDATE ON RANK %4i   MSG %4i   KEFF %16.12f \n",
//              Comm::rank(), msg, k);

  d_is_updated = true;
  if (soft_equiv(k, d_keff) && store_k)
  {
    // Do nothing, since these responses are already stored.
    d_is_updated = false;
  }
  else if (soft_equiv(k, d_keff_1) && store_k) // may want soft equiv
  {
    // Swap the stored eigenvalues
    d_keff_1 = d_keff;
    d_keff   = k;

    // Swap the responses
    std::swap(d_responses, d_responses_1);

    ++d_number_swapped_evaluations;
  }
  else
  {
    // Swap the stored eigenvalues
    d_keff_1 = d_keff;
    d_keff   = k;

    // Swap if we're storing
    if (store_k) std::swap(d_responses, d_responses_1);

    // Update sources for the new eigenvalue
    for (int i = 0; i < d_sources.size(); i++)
    {
      d_sources[i]->update(k);
      d_responses[i]->clear();
    }

    // Use the simple updater.
    update_explicit_work_share();

    ++d_number_unique_evaluations;
  }

  // Must go back to world, for which the global
  // roots now have the updated response.
  Comm::set(serment_comm::world);

  d_response_time += (Comm::wtime() - et);

  return d_is_updated;
}

//----------------------------------------------------------------------------//
bool ResponseServer::is_updated() const
{
  return d_is_updated;
}

//----------------------------------------------------------------------------//
ResponseServer::size_t ResponseServer::number_unique_evaluations() const
{
  return d_number_unique_evaluations;
}

//----------------------------------------------------------------------------//
ResponseServer::size_t ResponseServer::number_swapped_evaluations() const
{
  return d_number_swapped_evaluations;
}

//----------------------------------------------------------------------------//
double ResponseServer::response_time() const
{
  return d_response_time;
}

//----------------------------------------------------------------------------//
double ResponseServer::last_keff(const double k) const
{
  double val = k;
  if (soft_equiv(k, d_keff))
  {
    val = d_keff_1;
  }
  else if (soft_equiv(k, d_keff_1))
  {
    val = d_keff;
  }
  else
  {
    THROW("ERROR USING LAST KEFF.");
  }
  return val;
}

//----------------------------------------------------------------------------//
// IMPLEMENTATION
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
/*
 *  Simplest case of having each process do a predefined amount of
 *  work.  This first implementation assumes all responses are created
 *  equal.  Must be on local comm.
 *
 *  \todo This is the second use of the "split a vector among process"
 *        pattern; generalize it and put to Comm
 */
void ResponseServer::update_explicit_work_share()
{
  // Clear the responses.  This is needed since we use reduction to collect.
  // Ultimately, more sophisticated approached will be assessed.
  for (int i = 0; i < d_responses.size(); ++i)
    d_responses[i]->clear();

  size_t number_responses = 0;
  std::vector<size_t> number_per_process(Comm::size(), 0);

  // Local root
  // \todo use Comm::partition
  if (Comm::rank() == 0)
  {
    // Compute the number of unique local responses.
    number_responses = d_indexer->number_unique_moments();

    // Initial guess for responses per process and the remainder.
    int npp = number_responses / Comm::size();
    int remainder = number_responses - npp * Comm::size();

    // Assign the number per process, putting the extras on
    // the processes in reverse.  If there are more processes
    // than nodes, put the nodes on the first processes.
    int bound = Comm::size();
    if (!npp)
    {
      bound = number_responses;
      npp = 1;
      remainder = 0;
    }

    for (int i = 0; i < bound; i++)
      number_per_process[i] = npp;

    for (int i = 1; i <= remainder; i++)
      number_per_process[Comm::size() - i] += 1;

  }

  // Find my start and finish
  Comm::broadcast(&number_responses, 1, 0);
  Comm::broadcast(&number_per_process[0], number_per_process.size(), 0);
  size_t start = 0;
  for (int i = 0; i < Comm::rank(); i++)
    start += number_per_process[i];
  size_t finish = start + number_per_process[Comm::rank()];

  // Loop over all of my unique local moments
  for (size_t index_ul = start; index_ul < finish; index_ul++)
  {

    const ResponseIndex index_r =
      d_indexer->response_index_from_unique_local(index_ul);

    // Local unique node index
    size_t node_ul = d_nodes->unique_local_index_from_unique_global(index_r.node);

    // Compute responses
    Assert(node_ul < d_responses.size());
    Assert(node_ul < d_sources.size());

    d_sources[node_ul]->compute(d_responses[node_ul], index_r);

  }
  Comm::global_barrier();

  // A simple way to gather the results on 0 is to reduce on the
  // arrays of each nodal response.  Note, this probably makes the
  // best sense for small numbers of local processes.
  if (Comm::size() > 0)
  {
    for (size_t n = 0; n < d_sources.size(); n++)
    {
      int number_moments = d_responses[n]->size();
      int number_surfaces = d_responses[n]->number_surfaces();
      int number_pins = d_responses[n]->number_pins();
      for (size_t in = 0; in < number_moments; in++)
      {
        Comm::sum(&d_responses[n]->boundary_response(0, in), number_moments, 0);
        Comm::sum(&d_responses[n]->leakage_response(0, in), number_surfaces, 0);
        if (number_pins)
          Comm::sum(&d_responses[n]->pin_power(0, in), number_pins, 0);
      }
      Comm::sum(&d_responses[n]->fission_response(0), number_moments, 0);
      Comm::sum(&d_responses[n]->absorption_response(0), number_moments, 0);
      Comm::sum(&d_responses[n]->nodal_power(0), number_moments, 0);
    }
  }

}

} // end namespace erme_response

//----------------------------------------------------------------------------//
//              end of file ResponseServer.cc
//----------------------------------------------------------------------------//
