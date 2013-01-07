//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ResponseMatrix.cc
 *  @brief  ResponseMatrix member definitions
 *  @author Jeremy Roberts
 *  @date   Sep 2, 2012
 */
//---------------------------------------------------------------------------//

#include "ResponseMatrix.hh"

namespace erme
{

//---------------------------------------------------------------------------//
ResponseMatrix::ResponseMatrix(SP_nodelist nodes,
                               SP_indexer indexer,
                               SP_server server)
  : ResponseOperator(nodes, indexer, server)
  , Matrix(indexer->number_local_moments(),
           indexer->number_local_moments())
{

  /*
   *  Build number of nonzeros
   *
   *  In parallel, PETSc defines matrices in terms of local rows.  To
   *  preallocate a parallel matrix, the number of nonzeros per local
   *  row must be specified.  This is broken into "on diagonal" and
   *  "off diagonal" blocks.  The on diagonal blocks contain those
   *  elements that are within the local row/column block.  Off diagonal
   *  terms are within the local rows but are in columns that don't
   *  correspond to the local rows.  For the response matrix, there
   *  are no off diagonal terms.
   */
  vec_int nnz_on_diag(indexer->number_local_moments(),  0);
  vec_int nnz_off_diag(indexer->number_local_moments(), 0);

  int m = 0;
  for (size_t n = nodes->lower_bound(); n < nodes->upper_bound(); n++)
  {
    size_t un = nodes->unique_global_index_from_global(n);
    int size = indexer->number_node_moments(un);
    for (int i = 0; i < size; i++, m++)
      nnz_on_diag[m] = size;
  }

  // Preallocate.  This also computes the bounds and such.
  preallocate(nnz_on_diag, nnz_off_diag);

}

//---------------------------------------------------------------------------//
void ResponseMatrix::update()
{
  using std::cout;
  using std::endl;

  // Offset for a block.  Starts at this matrix's lower bound.
  int offset = lower_bound();

  // Loop through nodes
  for (int n = 0; n < d_nodes->number_local_nodes(); n++)
  {

    // Get response
    SP_response r = d_server->response(n);

    //r->display();

    // Block indices
    vec_int indices(r->size(), offset);
    for (size_type i = 0; i < indices.size(); i++)
      indices[i] += i;

    // Insert each column (corresponding to an incident moment)
    for (int in = 0; in < r->size(); in++)
    {

      int col = in + offset;
      insert_values(r->size(), &indices[0],
                    1, &col,
                    &r->boundary_response(0, in));
    }

    offset += indices.size();
  }

  // Assemble after finishing
  assemble();

}

} // end namespace erme

//---------------------------------------------------------------------------//
//              end of file ResponseMatrix.cc
//---------------------------------------------------------------------------//
