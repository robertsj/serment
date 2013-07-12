//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Jacobian.cc
 *  @brief Jacobian member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "Jacobian.hh"

namespace erme_solver
{

//----------------------------------------------------------------------------//
Jacobian::Jacobian(SP_server            server,
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
    d_MR = new OperatorMR(d_R, d_M);
  }
  Require(d_eps > 0.0);

  d_m = d_R->number_global_rows();
  d_m_full = d_m;
  if (serment_comm::Comm::is_last()) d_m_full += 2;

  // allocate etc.
  d_matrix = new Shell(d_m_full, *this);

  // working vector
  d_fd_MR = new Vector(d_m, 0.0);
}

//----------------------------------------------------------------------------//
void Jacobian::multiply(Vector &f, Vector &fp_times_f)
{
  Assert(f.global_size() == d_MR->number_global_rows() + 2);

  using serment_comm::Comm;

  // The Jacobian has the form
  //   | (M*R-lambda*I)   M*R_k*J              -J |
  //   | (F - k*L)        (F_k-k*L_k-L)*J       0 |
  //   | -J'              0                     0 |
  // which is computed in pieces. We approximate M*R_k*J as
  //     d/dk(M*R(k)*J) ~ M*(R(k+eps)-R(k))/eps * J
  // which is computed and stored upon updating the Jacobian. Similarly,
  //     d/dk((F - k*L)*f_J) = [F(k+eps)-F(k) - k*(L(k+eps)+L(k))]/eps * f_J
  //                         - L(k)*f_J
  // Note that "L" as written here the Jacobian refers to  a total "loss"
  // operator including absorption and global leakage.

  // The current moment vector
  Vector x_J(*d_x, d_m);

  // Residual components
  Vector f_J(f, d_m);
  double f_k = 0;
  double f_lambda = 0;

  // Outgoing vector components
  Vector fp_times_f_J(fp_times_f, d_m);
  double fp_times_f_k = 0.0;
  double fp_times_f_lambda = 0.0;

  // last process is assigned the extra unknowns, so broadcast to others
  if (Comm::is_last())
  {
    f_k      = f[d_m_full - 1];
    f_lambda = f[d_m_full - 2];
  }
  Comm::broadcast(&f_k,       Comm::last());
  Comm::broadcast(&f_lambda,  Comm::last());

  // Outgoing current component
  d_MR->multiply(f_J, fp_times_f_J);
  fp_times_f_J.add_a_times_x(-d_lambda, f_J);
  fp_times_f_J.add_a_times_x(-f_lambda, x_J);
  fp_times_f_J.add_a_times_x(f_k, *d_fd_MR);

  // Outgoing k component
  fp_times_f_k = d_F->dot(f_J) - d_k * (d_A->dot(f_J) + d_L->leakage(f_J))
               + d_fd_FAL * f_k;

  // Outgoing lambda component
  fp_times_f_lambda = -x_J.dot(f_J);

  if (Comm::is_last())
  {
    fp_times_f[d_m_full-2] = fp_times_f_k;
    fp_times_f[d_m_full-1] = fp_times_f_lambda;
  }

  return;
}

//----------------------------------------------------------------------------//
void Jacobian::multiply_transpose(Vector &v_in, Vector &v_out)
{
  THROW("NOT IMPLEMENTED");
}

//----------------------------------------------------------------------------//
void Jacobian::update(SP_vector x)
{
  using serment_comm::Comm;
  d_x = x;
  if (Comm::is_last())
  {
    d_k      = (*d_x)[d_m_full - 1];
    d_lambda = (*d_x)[d_m_full - 2];
  }
  Comm::broadcast(&d_k, Comm::last());
  Comm::broadcast(&d_lambda, Comm::last());

  // Insert unknown in current-sized vector
  Vector x_J(*d_x, d_m);

  // Compute the finite differenced components
  Vector fd_MR_tmp(d_m, 0.0);

  // for initial keff
  d_server->update(d_k);
  d_MR->multiply(x_J, fd_MR_tmp);
  double gain_1 = d_F->dot(x_J);
  double loss_1 = d_A->dot(x_J) + d_L->leakage(x_J);

  // for perturbed keff
  d_server->update(d_k + d_eps);
  d_MR->multiply(x_J, *d_fd_MR);
  d_fd_MR->add_a_times_x(-1.0, fd_MR_tmp);
  d_fd_MR->scale(1.0/d_eps);

  double gain_2 = d_F->dot(x_J);
  double loss_2 = d_A->dot(x_J) + d_L->leakage(x_J);
  // (d/dk)[ gain(k) - k*loss(k) ] ~ (g2-g1)/eps - k*(l2-l1/eps - l1
  d_fd_FAL = (gain_2-gain_1)/d_eps - d_k*(loss_2-loss_1)/d_eps - loss_1;

}

} // end namespace erme_solver

