//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   NonlinearResidual.i.hh
 *  @brief  NonlinearResidual.i
 *  @author Jeremy Roberts
 *  @date   Oct 1, 2012
 */
//---------------------------------------------------------------------------//

#ifndef erme_solver_NONLINEARRESIDUAL_I_HH_
#define erme_solver_NONLINEARRESIDUAL_I_HH_

#include "comm/Comm.hh"
#include <cmath>

namespace erme_solver
{

//---------------------------------------------------------------------------//
inline double NonlinearResidual::
compute_norm(Vector &x)
{
  // Break this into B, k, and l
  Vector B(d_R->number_local_columns(), 0.0);

  // k and lambda on last process
  double k = 0;
  double l = 0;
  if (serment_comm::Comm::rank() == serment_comm::Comm::size() - 1)
  {
    k = B[B.local_size()-2];
    l = B[B.local_size()-1];
  }

  // Broadcast k and lambda to all
  serment_comm::Comm::broadcast(&k, 1, serment_comm::Comm::size() - 1);
  serment_comm::Comm::broadcast(&l, 1, serment_comm::Comm::size() - 1);

  // Return the norm
  return compute_norm(B, k, l);
}

//---------------------------------------------------------------------------//
inline double NonlinearResidual::
compute_norm(Vector &B, const double k, const double l)
{
  // Temporary vector for residuals
  Vector temp(d_R->number_local_columns(), 0.0);
  Vector F_B(d_R->number_local_columns(), 0.0);
  double F_k = 0;
  double F_l = 0;

  // Compute M * R * B
  d_R->multiply(B,    temp);
  d_M->multiply(temp, F_B);

  // Compute M*R*B - Lambda*B
  F_B.add_a_times_x(-l, B);

  // Compute ||M*R*B - Lambda*B||^2
  double norm_B_sq = std::pow(F_B.norm(), 2);

  // Compute square of (F - k * L) * B
  double norm_K_sq = std::pow(d_F->dot(B) -
                              k * (d_A->dot(B) + d_L->leakage(B)), 2);

  std::cout << " F  = " << d_F->dot(B)
            << " A  = " << d_A->dot(B)
            << " L  = " << d_L->leakage(B) << std::endl;

  // Compute square of 0.5( 1 - B' * B)
  double norm_L_sq = std::pow(0.5 - 0.5 * std::pow(B.norm(), 2), 2);

  std::cout << " F_B= " << std::sqrt(norm_B_sq)
            << " F_K= " << std::sqrt(norm_K_sq)
            << " F_L= " << std::sqrt(norm_L_sq) << std::endl;

  // Return final L2 norm of residual
  return std::sqrt(norm_B_sq + norm_K_sq + norm_L_sq);
}

} // end namespace erme_solver

#endif // erme_solver_NONLINEARRESIDUAL_I_HH_

//---------------------------------------------------------------------------//
//              end of file NonlinearResidual.i.hh
//---------------------------------------------------------------------------//
