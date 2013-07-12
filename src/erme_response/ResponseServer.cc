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

#include <iostream>

namespace erme_response
{

//----------------------------------------------------------------------------//
ResponseServer::ResponseServer(SP_nodelist  nodes,
                               SP_indexer   indexer,
                               std::string  dbname,
                               size_t       dborder)
  : d_nodes(nodes)
  , d_indexer(indexer)
{
  // Preconditions
  Require(d_nodes);
  Require(d_indexer);

  d_sources.resize(d_nodes->number_unique_local_nodes());
  d_responses.resize(d_nodes->number_unique_local_nodes());

  // Build the response database if requested
  if (dbname != "")
  {
    d_rfdb = ResponseDatabase::Create(dbname, dborder);
  }

  /*
   *  Create the sources and responses.  Note, these are
   *  done for *unique* nodes.  Thus, repeated nodes are
   *  computed only once locally.
   */
  ResponseSourceFactory builder;
  for (size_t node_ul = 0; node_ul < d_sources.size(); node_ul++)
  {
    // Build the sources
    size_t node_ug =  d_nodes->unique_global_index_from_unique_local(node_ul);
    d_sources[node_ul] = builder.build(nodes->unique_node(node_ug), d_indexer);
    Assert(d_sources[node_ul]);

    // Build the nodal response containers
    d_responses[node_ul] =
        new NodeResponse(d_indexer->number_node_moments(node_ug),
                         d_nodes->unique_node(node_ug)->number_surfaces());
  }

}

//----------------------------------------------------------------------------//
void ResponseServer::update(const double keff)
{
  // Preconditions
  Require(serment_comm::communicator == serment_comm::world);

  using std::cout;
  using std::endl;

  typedef serment_comm::Comm Comm;

  // Switch to local communicator.
  Comm::set(serment_comm::local);

  // Broadcast new keff to local communicator.  The local root
  // must have it from the global solve.
  double k = keff;
  Comm::broadcast(&k, 1, 0);

  // Update sources
  for (int i = 0; i < d_sources.size(); i++)
  {
    d_sources[i]->update(k);
    d_responses[i]->clear();
  }

  // Use the simple updater.
  update_explicit_work_share();

  // Must go back to world, for which the global
  // roots now have the updated response.
  Comm::set(serment_comm::world);

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
  using std::cout;
  using std::endl;

  typedef serment_comm::Comm Comm;

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
      for (size_t in = 0; in < number_moments; in++)
      {
        Comm::sum(&d_responses[n]->boundary_response(0, in), number_moments, 0);
        Comm::sum(&d_responses[n]->leakage_response(0, in), number_surfaces, 0);
      }
      Comm::sum(&d_responses[n]->fission_response(0), number_moments, 0);
      Comm::sum(&d_responses[n]->absorption_response(0), number_moments, 0);
    }
  }

}

} // end namespace erme_response

//----------------------------------------------------------------------------//
//              end of file ResponseServer.cc
//----------------------------------------------------------------------------//
