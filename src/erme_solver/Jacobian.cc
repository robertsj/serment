//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Jacobian.cc
 *  @brief Jacobian member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "Jacobian.hh"
#include "GlobalSolverBase.hh"
#include "utilities/SoftEquivalence.hh"
#include <cstdio>

#define COUT(c) std::cout << c << std::endl;

namespace erme_solver
{

using serment_comm::Comm;

//----------------------------------------------------------------------------//
Jacobian::Jacobian(SP_server            server,
                   SP_indexer           indexer,
                   SP_responsecontainer responses,
									 const double         eps)
  : d_server(server)
  , d_indexer(indexer)
  , d_eps(eps)
  , d_use_previous_keff(eps <= 0.0)
  , d_k(0.0)
  , d_lambda(0.0)
  , d_fd_FAL(0.0)
  , d_time(0.0)
{
  Require(d_server);
  Require(d_indexer);

  if (Comm::is_global())
  {
    Require(responses);

    Comm::set(serment_comm::global);

    d_M = responses->M; Ensure(d_M);
    d_R = responses->R; Ensure(d_R);
    d_L = responses->L; Ensure(d_L);
    d_A = responses->A; Ensure(d_A);
    d_F = responses->F; Ensure(d_F);
    d_MR = new OperatorMR(d_R, d_M);

    d_m = d_R->number_global_rows();
    d_m_full = d_m;
    if (Comm::is_last())
      d_m_full += 2;

    d_matrix = new Shell(d_m_full, *this);
    d_fd_MR = new Vector(d_m, 0.0);

    Comm::set(serment_comm::world);
  }
}

//----------------------------------------------------------------------------//
void Jacobian::multiply(Vector &f, Vector &fp_times_f)
{
  Require(f.global_size() == d_MR->number_global_rows() + 2);
  //COUT(" JACOBIAN MULTIPLY....")
  Comm::tic();

  Comm::set(serment_comm::global);

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
    f_k      = f[d_m_full - 2];
    f_lambda = f[d_m_full - 1];
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
    fp_times_f[d_m_full - 2] = fp_times_f_k;
    fp_times_f[d_m_full - 1] = fp_times_f_lambda;
  }

  Comm::set(serment_comm::world);
  d_time += Comm::toc();
  //COUT(" JACOBIAN MULTIPLY....DONE")
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
  Require(serment_comm::communicator == serment_comm::world);

  Comm::tic();
  Comm::set(serment_comm::global);

  d_x = x;
  if (Comm::is_last())
  {
    d_k      = (*d_x)[d_m_full - 2];
    d_lambda = (*d_x)[d_m_full - 1];
  }
  Comm::broadcast(&d_k, Comm::last());
  Comm::broadcast(&d_lambda, Comm::last());

  Comm::set(serment_comm::world);
  update_response(d_k);
  Comm::set(serment_comm::global);

  // insert unknown in current-sized vector
  Vector x_J(*d_x, d_m);

  // compute the finite differenced components
  Vector fd_MR_tmp(d_m, 0.0);
  //   for initial keff
  d_MR->multiply(x_J, fd_MR_tmp);
  double gain_1 = d_F->dot(x_J);
  double loss_1 = d_A->dot(x_J) + d_L->leakage(x_J);
  //   for perturbed keff
  Comm::set(serment_comm::world);
  if (d_use_previous_keff)
  {
    d_eps = d_server->last_keff(d_k) - d_k;
    Assert(!detran_utilities::soft_equiv(d_eps, 0.0));
  }
  update_response(d_k + d_eps);
  Comm::set(serment_comm::global);
  d_MR->multiply(x_J, *d_fd_MR);
  double gain_2 = d_F->dot(x_J);
  double loss_2 = d_A->dot(x_J) + d_L->leakage(x_J);

  //   result
  d_fd_MR->add_a_times_x(-1.0, fd_MR_tmp);
  d_fd_MR->scale(1.0/d_eps);
  d_fd_FAL = (gain_2-gain_1)/d_eps - d_k*(loss_2-loss_1)/d_eps - loss_1;

  // return original responses
  Comm::set(serment_comm::world);
  update_response(d_k);

  d_time += Comm::toc();
}

//----------------------------------------------------------------------------//
void Jacobian::update_response(const double keff)
{
  Require(serment_comm::communicator == serment_comm::world);
  // alert the workers to update
  int msg = 1234; //GlobalSolverBase::CONTINUE;

  int ierr = serment_comm::Comm::broadcast(&msg, 1, 0);
  Assert(!ierr);

  // update the server and responses
  if (d_server->update(keff))
  {
    if (Comm::is_global())
    {
      d_R->update();
      d_F->update();
      d_A->update();
      d_L->update();
    }
  }
}

} // end namespace erme_solver

//----------------------------------------------------------------------------//
//              end of file Jacobian.cc
//----------------------------------------------------------------------------//

