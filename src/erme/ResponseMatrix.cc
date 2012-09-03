//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseMatrix.cc
 * \brief  ResponseMatrix member definitions
 * \author Jeremy Roberts
 * \date   Sep 2, 2012
 */
//---------------------------------------------------------------------------//

#include "ResponseMatrix.hh"

namespace erme
{

ResponseMatrix::ResponseMatrix(SP_nodelist nodes,
                               SP_indexer indexer,
                               SP_server server)
  : ResponseOperator(nodes, indexer, server)
  , linear_algebra::Matrix(indexer->number_local_moments(),
                           indexer->number_local_moments(),
                           // We have at most a "block size" per row, all
                           // within the local range.
                           vec_int(indexer->number_local_moments(),
                                   indexer->number_local_moments()),
                           // There are no off block-diagonal entries.
                           vec_int(indexer->number_local_moments(),
                                   0))
{

}

void ResponseMatrix::update(const double keff)
{

  // Update the responses
  d_server->update(keff);

  // Loop through nodes
  for (int n = 0; n < d_nodes->number_local_nodes(); n++)
  {
    // Get response
    SP_response r = d_server->response(n);

    // Block indices
    vec_int indices(r->size(), lower_bound());
    for (size_type i = 0; i < indices.size(); i++)
      indices[i] += i;

    // Insert each column (corresponding to an incident moment)
    for (int in = 0; in < r->size(); in++)
    {
      insert_values(r->size(), &indices[0], 1, &in,
                    &r->boundary_response(0, in));
    }
  }

  // Assemble after finishing
  assemble();

}

} // end namespace erme

//---------------------------------------------------------------------------//
//              end of file ResponseMatrix.cc
//---------------------------------------------------------------------------//
