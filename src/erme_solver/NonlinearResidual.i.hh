//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  NonlinearResidual.i.hh
 *  @brief NonlinearResidual inline member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef erme_solver_NONLINEARRESIDUAL_I_HH_
#define erme_solver_NONLINEARRESIDUAL_I_HH_

#include "comm/Comm.hh"
#include <cmath>

namespace erme_solver
{

//----------------------------------------------------------------------------//
inline void NonlinearResidual::
evaluate(Vector &x, Vector &f)
{
  using serment_comm::Comm;

  // Create temporary of smaller size with f's memory.
  Vector x_J(x, d_MR->number_local_rows());
  Vector f_J(f, d_MR->number_local_rows());

  // Broadcast eigenvalues k and lambda from last process to all
  double k = 0.0;
  double l = 0.0;
  int    m = x.local_size();
  if (Comm::is_last())
  {
    k = x[m-2];
    l = x[m-1];
  }
  Comm::broadcast(&k, 1, Comm::last());
  Comm::broadcast(&l, 1, Comm::last());

//  {
//    d_MR->display(d_MR->BINARY, "MR.out");
//    d_L->leakage_vector().display(d_MR->BINARY, "L.out");
//    d_F->display(d_MR->BINARY, "F.out");
//    d_A->display(d_MR->BINARY, "A.out");
//    x.display(d_MR->BINARY, "X.out");
//  }

  // Compute MR*J - l*J
  d_MR->multiply(x_J, f_J);
  f_J.add_a_times_x(-l, x_J);

  // Compute f_k = gains - k*losses
  double f_k = d_F->dot(x_J) - k * (d_A->dot(x_J) + d_L->leakage(x_J));

  // Compute f_l
  double f_l = 0.5 - 0.5 * std::pow(x_J.norm(), 2);

  if (serment_comm::Comm::is_last())
  {
    f[m-2] = f_k;
    f[m-1] = f_l;
  }
}

//----------------------------------------------------------------------------//
inline double NonlinearResidual::
compute_norm(Vector &x)
{
  Vector f(x.local_size(), 0.0);
  evaluate(x, f);
  return f.norm(f.L2);
}

//----------------------------------------------------------------------------//
inline double NonlinearResidual::
compute_norm(Vector &J, const double k, const double l)
{
  size_t m = J.local_size();
  if (serment_comm::Comm::is_last())
    m += 2;
  Vector x(m, 0);
  for (size_t i = 0; i < J.local_size(); ++i)
    x[i] = J[i];
  if (serment_comm::Comm::is_last())
  {
    x[m-2] = k;
    x[m-1] = l;
  }
  return compute_norm(x);
}



} // end namespace erme_solver

#endif // erme_solver_NONLINEARRESIDUAL_I_HH_

//----------------------------------------------------------------------------//
//              end of file NonlinearResidual.i.hh
//----------------------------------------------------------------------------//
