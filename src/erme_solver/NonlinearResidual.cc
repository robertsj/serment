//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  NonlinearResidual.cc
 *  @brief NonlinearResidual
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "NonlinearResidual.hh"
#include "utilities/DBC.hh"

namespace erme_solver
{

using serment_comm::Comm;

//----------------------------------------------------------------------------//
NonlinearResidual::NonlinearResidual(SP_server server, SP_responses responses)
  : d_server(server)
  , d_responses(responses)
{
  Require(d_server);
  if (serment_comm::Comm::is_global())
  {
    Require(d_responses);
    d_MR = new OperatorMR(d_responses->R, d_responses->M);
  }
}

//----------------------------------------------------------------------------//
void NonlinearResidual::evaluate(Vector *x, Vector *f)
{
  Require(serment_comm::communicator == serment_comm::world);

  /*
   *  The structure of this call allows for it to be called simultaneously
   *  from the the global communicator or from within a global subset.  For
   *  the latter case, the non global block in which a single update is
   *  called is *not* called, and the client must call this somewhere else.
   */

  if (Comm::is_global())
  {
    Require(x);
    Require(f);

    // Switch to global communicator to handle unknowns
    Comm::set(serment_comm::global);

    // Create temporary of smaller size with f's memory.
    Vector x_J(*x, d_MR->number_local_rows());
    Vector f_J(*f, d_MR->number_local_rows());

    // Broadcast eigenvalues k and lambda from last process to all
    double k = 0.0;
    double l = 0.0;
    int    m = x->local_size();
    if (Comm::is_last())
    {
      k = (*x)[m-2];
      l = (*x)[m-1];
    }
    Comm::broadcast(&k, 1, Comm::last());
    Comm::broadcast(&l, 1, Comm::last());

    // Update server on the world communicator
    Comm::set(serment_comm::world);
    d_server->update(k);
    Comm::set(serment_comm::global);

    // Update the responses with the new eigenvalue
    d_responses->update();

    // boundary component: residual of the current eigenvalue problem
    d_MR->multiply(x_J, f_J);
    f_J.add_a_times_x(-l, x_J);

    // keff component: gains - k * losses
    double f_k = d_responses->F->dot(x_J) -
                 k * (d_responses->A->dot(x_J) + d_responses->L->leakage(x_J));

    // lambda component: ensures normalization of boundary unknowns
    double f_l = 0.5 - 0.5 * std::pow(x_J.norm(), 2);

    if (Comm::is_last())
    {
      (*f)[m-2] = f_k;
      (*f)[m-1] = f_l;
    }

    // Reset to world
    Comm::set(serment_comm::world);
  }
  else
  {
    // Update server.  Note that keff is broadcasted within, so dummy is fine.
    d_server->update(0.0);
  }
}

//----------------------------------------------------------------------------//
double NonlinearResidual::compute_norm(Vector *x)
{
  Require(serment_comm::communicator == serment_comm::world);
  if (Comm::is_global()) Require(x);

  SP_vector f;
  if (Comm::is_global())
    f = new Vector(x->local_size(), 0.0);

  evaluate(x, f.bp());

  double norm_f = 0.0;
  if (Comm::is_global())
    norm_f = f->norm(f->L2);
  Comm::broadcast(&norm_f, 1, 0);

  return norm_f;
}

//----------------------------------------------------------------------------//
double NonlinearResidual::compute_norm(Vector       *J,
                                       const double  k,
                                       const double  l)
{
  Require(serment_comm::communicator == serment_comm::world);
  if (Comm::is_global()) Require(J);

  SP_vector x;

  if (Comm::is_global())
  {
    Comm::set(serment_comm::global);
    size_t m = J->local_size();
    if (Comm::is_last())
      m += 2;
    x = new Vector(m, 0.0);
    for (size_t i = 0; i < J->local_size(); ++i)
      (*x)[i] = (*J)[i];
    if (Comm::is_last())
    {
      (*x)[m-2] = k;
      (*x)[m-1] = l;
    }
    Comm::set(serment_comm::world);
  }

  return compute_norm(x.bp());
}

} // end namespace erme_solver

//----------------------------------------------------------------------------//
//              end of file NonlinearResidual.cc
//----------------------------------------------------------------------------//



