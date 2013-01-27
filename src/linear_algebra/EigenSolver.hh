//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   EigenSolver.hh
 *  @author robertsj
 *  @date   Sep 5, 2012
 *  @brief  EigenSolver class definition.
 */
//---------------------------------------------------------------------------//

#ifndef linear_algebra_EIGENSOLVER_HH_
#define linear_algebra_EIGENSOLVER_HH_EIGENSOLVER_HH_

#include "Vector.hh"
#include "MatrixBase.hh"
#include "comm/Comm.hh"
#include "utilities/DBC.hh"
#include "utilities/SP.hh"
#include "slepc.h"
#include <vector>
#include <string>

namespace linear_algebra
{

class EigenSolver
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<EigenSolver>     SP_solver;
  typedef MatrixBase::SP_matrix                 SP_matrix;
  typedef Vector::SP_vector                     SP_vector;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param A  Pointer to left hand operator
   *  @param P  Pointer to right hand operator (possibly null)
   */
  EigenSolver(SP_matrix A,
              SP_matrix B = SP_matrix(0));

  /**
   *  @brief Solve the eigenvalue problem
   *  @param x  Eigenvector
   *  @return   Dominant eigenvalue
   */
  double solve(SP_vector x);

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// SLEPc eigensolver
  EPS d_solver;
  /// Left side
  SP_matrix d_A;
  /// Right side
  SP_matrix d_B;
  /// Temporary vector for imaginary component
  SP_vector d_x_imaginary;
  /// Eigenvalue (real part)
  double d_lambda;
  /// Eigenvalue (imaginary part)
  double d_lambda_imaginary;
  /// Number of iterations
  int d_number_iterations;
  /// Maximum number of iterations
  int d_maximum_iterations;
  /// Tolerance
  double d_tolerance;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

};

} // end namespace linear_algebra

#endif /* linear_algebra_EIGENSOLVER_HH_EIGENSOLVER_HH_ */

//---------------------------------------------------------------------------//
//              end of file EigenSolver.hh
//---------------------------------------------------------------------------//
