//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ResponseServer.cc
 *  @author robertsj
 *  @date   Aug 31, 2012
 *  @brief  ResponseServer class definition.
 */
//---------------------------------------------------------------------------//

#include "ResponseServer.hh"
#include "ResponseSourceFactory.hh"
#include "comm/Comm.hh"

#include <iostream>

namespace erme_response
{

ResponseServer::ResponseServer(SP_nodelist  nodes,
                               SP_indexer   indexer,
                               std::string  dbname)
  : d_nodes(nodes)
  , d_indexer(indexer)
  , d_sources(nodes->number_local_nodes())
  , d_responses(nodes->number_local_nodes())
{
  // Preconditions
  Require(nodes);
  Require(indexer);

  // Build the response database if requested
  if (dbname != "")
  {
    d_rfdb = ResponseDatabase::Create(dbname);
  }

  ResponseSourceFactory builder;
  for (size_t n = 0; n < d_sources.size(); n++)
  {
    // Build the sources
    size_t n_global = nodes->global_index(n);
    d_sources[n] = builder.build(nodes->node(n_global));
    Ensure(d_sources[n]);

    // Build the nodal response containers
    d_responses[n] =
        new NodeResponse(d_indexer->number_node_moments(n_global),
                         d_nodes->node(n_global)->number_surfaces());
    Ensure(d_responses[n]);
  }

}

/*
 * Currently, I expend the updates to be called from the global
 * communicator, since they want
 *
 */
void ResponseServer::update(const double keff)
{
  // Preconditions
  Require(serment_comm::communicator == serment_comm::world);

  using std::cout;
  using std::endl;

  typedef serment_comm::Comm Comm;

  // Switch to local communicator.
  Comm::set(serment_comm::local);
  double k = keff;

  // Broadcast new keff to local communicator.  The local root
  // must have it from the global solve.
  Comm::broadcast(&k, 1, 0);

  // Update sources
  for (int i = 0; i < d_sources.size(); i++)
    d_sources[i]->update(k);

  // Use the simple updater.
  update_explicit_work_share();

  // Must go back to world, for which the global
  // roots now have the updated response.
  Comm::set(serment_comm::world);

}

//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//

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

  // total number
  size_t number_responses = 0;
  std::vector<size_t> number_per_process(Comm::size(), 0);

  // Local root \todo use Comm::partition
  if (Comm::rank() == 0)
  {
    number_responses = d_indexer->number_local_moments();

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

 // Broadcast the number of nodes in the problem
 Comm::broadcast(&number_responses, 1, 0);

 // Broadcast the number of nodes per process.
 Comm::broadcast(&number_per_process[0], number_per_process.size(), 0);

 // Find my start and finish
 size_t start = 0;
 for (int i = 0; i < Comm::rank(); i++)
   start += number_per_process[i];
 size_t finish = start + number_per_process[Comm::rank()];

 // Loop over all of my local moments
 for (size_t index_l = start; index_l < finish; index_l++)
 {

   const ResponseIndex index_r = d_indexer->response_index(index_l);

   // Local node index
   int node_l = d_nodes->local_index(index_r.node);

   // Compute responses
   Assert(node_l < d_responses.size());
   Assert(node_l < d_sources.size());

   d_sources[node_l]->compute(d_responses[node_l], index_r);

 }
 Comm::global_barrier();

 // A simple way to gather the results on 0 is to reduce on the
 // arrays of each nodal response.  Note, this probably makes the
 // best sense for small numbers of local processes
 for (size_t n = 0; n < d_sources.size(); n++)
 {
   int number_moments = d_responses[n]->size();
   int number_surfaces = d_responses[n]->number_surfaces();
   for (size_t in = 0; in < number_moments; in++)
   {
//     if (Comm::world_rank() == 0)
//     {
//       cout << " was " << d_responses[n]->boundary_response(0, in) << endl;
//     }

     Comm::sum(&d_responses[n]->boundary_response(0, in), number_moments,  0);
     Comm::sum(&d_responses[n]->leakage_response(0, in),  number_surfaces, 0);

//     if (Comm::world_rank() == 0)
//     {
//       cout << " now " << d_responses[n]->boundary_response(0, in) << endl;
//     }
   }
   Comm::sum(&d_responses[n]->fission_response(0),    number_moments, 0);
   Comm::sum(&d_responses[n]->absorption_response(0), number_moments, 0);
 }

}

} // end namespace erme_response


