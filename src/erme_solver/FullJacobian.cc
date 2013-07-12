//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  FullJacobian.cc
 *  @brief FullJacobian member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "FullJacobian.hh"
#include "utilities/MathUtilities.hh"

namespace erme_solver
{

//----------------------------------------------------------------------------//
FullJacobian::FullJacobian(SP_server            server,
                           SP_responsecontainer responses,
                           const double         eps)
  : d_server(server)
  , d_eps(eps)
  , d_k(0.0)
  , d_lambda(0.0)
  , d_fd_FAL(0.0)
{
  Require(d_server);
  if (serment_comm::Comm::is_global())
  {
    Require(responses);
    d_M = responses->M; Ensure(d_M);
    d_R = responses->R; Ensure(d_R);
    d_L = responses->L; Ensure(d_L);
    d_A = responses->A; Ensure(d_A);
    d_F = responses->F; Ensure(d_F);
  }
  Require(d_eps > 0.0);

  d_m = d_R->number_global_rows();
  d_m_full = d_m;
  if (serment_comm::Comm::is_last()) d_m_full += 2;

  /*
   *  Number of the nonzeros should be the number of moments of all local
   *  nodes plus one for the k-column plus one for the lambda-column.  For
   *  the last process, the k-row is full except the last element, and
   *  the lambda-row is full except for the last two elements
   *
   */
  vec_int nnz(d_m_full,    0);
  vec_int nnz_od(d_m_full, 0);

  // allocate etc.
  d_matrix = new linear_algebra::Matrix(d_m_full, d_m_full, nnz, nnz_od);

  // working vector
  d_fd_MR = new Vector(d_m, 0.0);
}


//----------------------------------------------------------------------------//
void FullJacobian::update(SP_vector x)
{
  using serment_comm::Comm;
  using detran_utilities::range;

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
  vec_int J_idx_v  = range<int>(d_R->lower_bound(), d_R->upper_bound(), false);
  int* J_idx       = &J_idx_v[0];
  int k_idx[]      = {d_R->number_global_columns()};
  int lambda_idx[] = {d_R->number_global_columns()};

  //--------------------------------------------------------------------------//
  // (M*R-LAMBDA*I) BLOCK
  //--------------------------------------------------------------------------//

  // create a new matrix with the product M*R
  linear_algebra::Matrix::SP_matrix MR = d_M->multiply(d_R);
  // copy MR to the Jacobian matrix
  d_matrix->insert_values(MR);
  // insert -lambda
  d_matrix->assemble(d_matrix->FLUSH);
  for (int i = d_R->lower_bound(); d_R->upper_bound(); ++i)
  {
    int diag[] = {i};
    double val[] = {-d_lambda};
    d_matrix->insert_values(1, diag, 1, diag, val, d_matrix->ADD);
  }
  d_matrix->assemble(d_matrix->FLUSH);

  //--------------------------------------------------------------------------//
  // K-ROW (EXCLUDING K-COLUMN)
  //--------------------------------------------------------------------------//

  Vector k_row(d_m, 0.0);
  k_row.add(*d_F);
  k_row.add_a_times_x(-d_k, *d_A);
  k_row.add_a_times_x(-d_k, d_L->leakage_vector());
  Vector::SP_vector k_row_seq = k_row.collect_on_root(Comm::last());
  if (Comm::rank() == Comm::last())
    d_matrix->insert_values(1, k_idx, d_m, J_idx, &(*k_row_seq)[0]);

  //--------------------------------------------------------------------------//
  // FINITE DIFFERENCE (I.E. K) COLUMN
  //--------------------------------------------------------------------------//

  // temporary for (d/dk)MR(k)J
  Vector fd_MR_tmp(d_m, 0.0);

  // initial keff
  d_server->update(d_k);
  MR->multiply(x_J, fd_MR_tmp);
  double gain_1 = d_F->dot(x_J);
  double loss_1 = d_A->dot(x_J) + d_L->leakage(x_J);

  // perturbed keff
  d_server->update(d_k + d_eps);
  MR->multiply(x_J, *d_fd_MR);

  // compute and insert perturbed column
  d_fd_MR->add_a_times_x(-1.0, fd_MR_tmp);
  d_fd_MR->scale(1.0/d_eps);
  double gain_2 = d_F->dot(x_J);
  double loss_2 = d_A->dot(x_J) + d_L->leakage(x_J);
  d_fd_FAL = (gain_2-gain_1)/d_eps - d_k*(loss_2-loss_1)/d_eps - loss_1;
  d_matrix->insert_values(d_m, J_idx, 1, k_idx, &(*d_fd_MR)[0]);
  if (Comm::is_last())
  {
    double vals[] = {d_fd_FAL};
    d_matrix->insert_values(1, k_idx, 1, k_idx, vals);
  }

  //--------------------------------------------------------------------------//
  // LAST TWO ROWS AND LAST COLUMN
  //--------------------------------------------------------------------------//

  // negate the current and collect on last process
  x_J.scale(-1.0);
  Vector::SP_vector lambda_row_seq = x_J.collect_on_root(Comm::last());

  // insert -J in last column and row
  d_matrix->insert_values(d_m, J_idx, 1, lambda_idx, &x_J[0]);
  if (Comm::rank() == Comm::last())
    d_matrix->insert_values(1, lambda_idx, d_m, J_idx, &x_J[0]);

}

} // end namespace erme_solver

