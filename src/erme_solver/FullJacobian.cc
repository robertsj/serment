//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  FullJacobian.cc
 *  @brief FullJacobian member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "FullJacobian.hh"
#include "GlobalSolverBase.hh"
#include "utilities/MathUtilities.hh"
#include <cstdio>

#define COUT(c) std::cout << c << std::endl;

namespace erme_solver
{

using serment_comm::Comm;
using detran_utilities::range;

//----------------------------------------------------------------------------//
FullJacobian::FullJacobian(SP_server            server,
                           SP_indexer           indexer,
                           SP_responsecontainer responses,
                           const double         eps,
                           bool                 fill_diag,
                           bool                 include_fd)
  : d_server(server)
  , d_indexer(indexer)
  , d_eps(eps)
  , d_fill_diag(fill_diag)
  , d_include_fd(include_fd)
  , d_k(0.0)
  , d_lambda(0.0)
  , d_fd_FAL(0.0)
{
  Require(d_server);
  Require(d_indexer);
  Require(d_eps > 0.0);

  if (serment_comm::Comm::is_global())
  {
    Require(responses);
    d_M = responses->M; Ensure(d_M);
    d_R = responses->R; Ensure(d_R);
    d_L = responses->L; Ensure(d_L);
    d_A = responses->A; Ensure(d_A);
    d_F = responses->F; Ensure(d_F);
    d_MR = new OperatorMR(d_R, d_M);
  }

  d_m = d_R->number_global_rows();
  d_m_full = d_m;
  if (serment_comm::Comm::is_last())
    d_m_full += 2;

  if (serment_comm::Comm::is_last())
  {
    std::cout << " I AM LAST " << std::endl;
  }
  else
  {
    THROW("why am I not last?");
  }

  d_matrix = new linear_algebra::Matrix(d_m_full, d_m_full);
  linear_algebra::Matrix &M =
      *dynamic_cast<linear_algebra::Matrix*>(d_matrix.bp());
  std::cout << " J size = " << d_matrix->number_global_columns()
            << " R size = " << d_R->number_global_columns() << std::endl;
  /*
   *  Number of the nonzeros should be the number of moments of all local
   *  nodes plus one for the k-column plus one for the lambda-column.  For
   *  the last process, the k-row is full except the last element, and
   *  the lambda-row is full except for the last two elements
   *
   */
  vec_int nnz(d_m_full,    0);
  vec_int nnz_od(d_m_full, 0);
  size_t i = 0;
  for (size_t gn = 0; gn < d_indexer->nodes()->number_local_nodes(); ++gn)
  {
    size_t ugn = d_indexer->nodes()->unique_global_index_from_global(gn);
    for (size_t m = 0; m < d_indexer->number_node_moments(ugn); ++m, ++i)
    {
      Assert(i < d_m_full);
      nnz[i]    = d_indexer->number_node_moments(ugn);
      nnz_od[i] = d_indexer->number_node_moments(ugn);
    }
  }
  if (serment_comm::Comm::is_last())
  {
    nnz[i  ] = d_m_full - 1;
    nnz[i+1] = d_m_full - 2;
    nnz_od[i]   = d_R->number_global_rows() + 2 - d_m_full;
    nnz_od[i+1] = d_R->number_global_rows() + 2 - d_m_full;
  }

  // allocate etc.
  M.preallocate(nnz, nnz_od);

  // working vector
  d_fd_MR = new Vector(d_m, 0.0);

}

//----------------------------------------------------------------------------//
void FullJacobian::multiply(Vector &v_in, Vector &v_out)
{
  v_out.copy(v_in);
  //d_matrix->multiply(v_in, v_out);
}

//----------------------------------------------------------------------------//
void FullJacobian::multiply_transpose(Vector &v_in, Vector &v_out)
{
  d_matrix->multiply_transpose(v_in, v_out);
}

//----------------------------------------------------------------------------//
void FullJacobian::update(SP_vector x)
{
  Require(serment_comm::communicator == serment_comm::world);

  // debug as identity
  Comm::set(serment_comm::global);
//  for (int i = d_matrix->lower_bound(); i < d_matrix->upper_bound(); ++i)
//  {
//    double val = 1.0;
//    d_matrix->insert_values(1, &i, 1, &i, &val);
//  }
//  d_matrix->assemble();
  //

  d_x = x;
  if (Comm::is_last())
  {
    d_k      = (*d_x)[d_m_full - 1];
    d_lambda = (*d_x)[d_m_full - 2];
  }
  Comm::broadcast(&d_k, Comm::last());
  Comm::broadcast(&d_lambda, Comm::last());
  Vector x_J(*d_x, d_m);

  // indices for the J-related, k-related, and lambda-related rows/columns
  vec_int J_idx_v      = range<int>(d_R->lower_bound(), d_R->upper_bound(), false);
  vec_int full_J_idx_v = range<int>(0, d_R->number_global_rows(), false);
  int* J_idx       = &J_idx_v[0];
  int* full_J_idx  = &full_J_idx_v[0];
  int k_idx[]      = {d_R->number_global_columns()  };
  int lambda_idx[] = {d_R->number_global_columns()+1};

  //--------------------------------------------------------------------------//
  // (M*R-LAMBDA*I) BLOCK
  //--------------------------------------------------------------------------//

  // create a new matrix with the product M*R
  linear_algebra::Matrix::SP_matrix MR = d_M->multiply(d_R);

  // copy MR to the Jacobian matrix
  d_matrix->insert_values(MR);
  // insert -lambda
  d_matrix->assemble(d_matrix->FLUSH);

  for (int i = d_R->lower_bound(); i < d_R->upper_bound(); ++i)
  {
    int diag[] = {i};
    double val[] = {-d_lambda};
    d_matrix->insert_values(1, diag, 1, diag, val, d_matrix->ADD);
  }
  d_matrix->assemble(d_matrix->FLUSH);

  //--------------------------------------------------------------------------//
  // K-ROW (EXCLUDING K-COLUMN)
  //--------------------------------------------------------------------------//

  Vector k_row(*d_A);
  k_row.add(d_L->leakage_vector());
  k_row.scale(-d_k);
  k_row.add(*d_F);
  Vector::SP_vector k_row_seq = k_row.collect_on_root(Comm::last());
  if (Comm::rank() == Comm::last())
  {
    d_matrix->insert_values(1, k_idx, d_m, full_J_idx, &(*k_row_seq)[0]);
  }

  //--------------------------------------------------------------------------//
  // FINITE DIFFERENCE (I.E. K) COLUMN
  //--------------------------------------------------------------------------//

  if (d_include_fd)
  {

    // temporary for (d/dk)MR(k)J
    Vector fd_MR_tmp(d_m, 0.0);
    //   initial keff
    Comm::set(serment_comm::world);
    d_server->update(d_k); // just the server, since this keff was alread done
    Comm::set(serment_comm::global);
    d_MR->multiply(x_J, fd_MR_tmp);
    double gain1 = d_F->dot(x_J);
    double loss1 = d_A->dot(x_J) + d_L->leakage(x_J);

    //   perturbed keff
    Comm::set(serment_comm::world);
    update_response(d_k + d_eps);
    Comm::set(serment_comm::global);
    d_MR->multiply(x_J, *d_fd_MR);
    double gain2 = d_F->dot(x_J);
    double loss2 = d_A->dot(x_J) + d_L->leakage(x_J);

    d_matrix->assemble();
    Comm::set(serment_comm::world);
    update_response(d_k);
    return;

    //   result
    d_fd_MR->add_a_times_x(-1.0, fd_MR_tmp);
    d_fd_MR->scale(1.0 / d_eps);
    d_fd_FAL = (gain2 - gain1) / d_eps - d_k * (loss2 - loss1) / d_eps - loss1;

    d_matrix->insert_values(d_m, J_idx, 1, k_idx, &(*d_fd_MR)[0]);
    if (Comm::is_last())
    {
      double vals[] = {d_fd_FAL};
      d_matrix->insert_values(1, k_idx, 1, k_idx, vals);
    }

  }
  else
  {
    if (d_fill_diag)
    {
      double diag_val = 1.0;
      d_matrix->insert_values(1, k_idx, 1, k_idx, &diag_val);
    }
  }

  //--------------------------------------------------------------------------//
  // LAST TWO ROWS AND LAST COLUMN
  //--------------------------------------------------------------------------//

  // negate the current and collect on last process
  x_J.scale(-1.0);
  Vector::SP_vector lambda_row_seq = x_J.collect_on_root(Comm::last());
  // insert -J in last column
  d_matrix->insert_values(d_m, J_idx, 1, lambda_idx, &x_J[0]);
  // insert -J in last row
  if (Comm::rank() == Comm::last())
  {
    d_matrix->insert_values(1, lambda_idx, d_m, full_J_idx,
                            &((*lambda_row_seq)[0]));
    if (d_fill_diag)
    {
      double diag_val = 1.0;
      d_matrix->insert_values(1, lambda_idx, 1, lambda_idx, &diag_val);
    }
  }
  // reset the sign, since this is our active unknown
  x_J.scale(-1.0);

  d_matrix->assemble();

  // return original responses
  Comm::set(serment_comm::world);
  update_response(d_k);
}

//----------------------------------------------------------------------------//
void FullJacobian::update_response(const double keff)
{
  Require(serment_comm::communicator == serment_comm::world);
  // alert the workers to update
  int msg = GlobalSolverBase::CONTINUE;
  serment_comm::Comm::broadcast(&msg, 1, 0);
  // update the server and responses
  d_server->update(keff);
  if (Comm::is_global())
  {
    d_R->update();
    d_F->update();
    d_A->update();
    d_L->update();
  }
}

} // end namespace erme_solver

//----------------------------------------------------------------------------//
//              end of file FullJacobian.cc
//----------------------------------------------------------------------------//
