//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  LeakageOperator.cc
 *  @brief LeakageOperator member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "LeakageOperator.hh"

namespace erme
{

//----------------------------------------------------------------------------//
LeakageOperator::LeakageOperator(SP_nodelist nodes,
                                 SP_indexer  indexer,
                                 SP_server   server)
  : ResponseOperator(nodes, indexer, server)
  , Matrix(nodes->number_local_surfaces(),
           indexer->number_local_moments())
  , d_global_leakage(nodes->number_local_surfaces(), 0.0)
  , d_L_times_moments(nodes->number_local_surfaces(), 0.0)
  , d_leakage_vector(indexer->number_local_moments(), 0.0)
{
  vec_int nnz_on_diag(d_nodes->number_local_surfaces(),  0);
  vec_int nnz_off_diag(d_nodes->number_local_surfaces(), 0);

  /*
   *  Note that the L operator has the form for a two 1D nodes
   *  with two moments per side:
   *
   *     --- --- --- --- --- --- --- ---
   *    | x | x | x | x | o | o | o | o |
   *     --- --- --- --- --- --- --- ---
   *    | x | x | x | x | o | o | o | o |
   *     --- --- --- --- --- --- --- ---
   *    | o | o | o | o | x | x | x | x |
   *     --- --- --- --- --- --- --- ---
   *    | o | o | o | o | x | x | x | x |
   *     --- --- --- --- --- --- --- ---
   *
   *  Per PETSc's matrix definitions, each row has just
   *  one nonzero on its diagonal.  The number of off diagonal
   *  nonzeros is the number of nodal moments minus one.
   *
   */

  int index_s = 0;
  int index_m = 0;

  // Loop over all nodes
  for (int n = d_nodes->lower_bound(); n < d_nodes->upper_bound(); n++)
  {
    // Unique node
    size_t un = nodes->unique_global_index_from_global(n);

    // Total moments for this node
    int size = d_indexer->number_node_moments(un);
    std::cout << " size = " << size << std::endl;

    // Loop over this nodes surfaces
    for (int s = 0; s < d_nodes->node(n)->number_surfaces(); ++s, ++index_s)
    {

      Assert(index_s < nnz_on_diag.size());
      nnz_on_diag[index_s] = 1;
      nnz_off_diag[index_s] = size - 1;// + 1;

      // Check for leakage
      if (d_nodes->neighbor(n, s).neighbor() == erme_geometry::Node::VACUUM)
      {
        Assert(index_s < d_global_leakage.local_size());
        d_global_leakage[index_s] = 1.0;
      }

    }

  }

  // Preallocate.  This also computes the bounds and such.
  preallocate(nnz_on_diag, nnz_off_diag);
  assemble();
}

//----------------------------------------------------------------------------//
void LeakageOperator::update()
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

    // Block indices for the rows.  These index into the
    // nodal surfaces.
    vec_int indices(r->number_surfaces(), offset);
    for (size_type i = 0; i < r->number_surfaces(); i++)
    {
      indices[i] += i;
    }

    // Converting the first nodal moment to the global
    // moment yields the offset for the column indices
    int ng = d_nodes->global_index_from_local(n);
    int column_offset = d_indexer->nodal_index_to_global(ng, 0);

    // Insert each column (corresponding to an incident moment)
    for (int in = 0; in < r->size(); in++)
    {
      int col = in + column_offset;
      insert_values(r->number_surfaces(), &indices[0],
                    1, &col,
                    &r->leakage_response(0, in));
    }

    offset += indices.size();
  }

  // Assemble after finishing
  assemble();

  d_L_times_moments.assemble();
  d_global_leakage.assemble();
}

//----------------------------------------------------------------------------//
double LeakageOperator::leakage(linear_algebra::Vector &x)
{
  Require(x.local_size() == number_local_columns());

  // Compute L * x
  multiply(x, d_L_times_moments);

  // Compute net global leakage  gL'*(L*x) = [1xs][s*n]*[n*s]
  // L'*gl = [n*s][s*1]
  double val = d_L_times_moments.dot(d_global_leakage);

  // \todo Could put warning for negative leakage
  return val;
}

//----------------------------------------------------------------------------//
const linear_algebra::Vector& LeakageOperator::leakage_vector()
{
  multiply_transpose(d_global_leakage, d_leakage_vector);
  return d_leakage_vector;
}

//----------------------------------------------------------------------------//
void LeakageOperator::display_leakage()
{
  std::cout << " GLOBAL LEAKAGE: " << std::endl;
  d_global_leakage.display();
}

} // end namespace erme

//----------------------------------------------------------------------------//
//              end of file LeakageOperator.cc
//----------------------------------------------------------------------------//
