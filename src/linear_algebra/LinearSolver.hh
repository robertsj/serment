//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   LinearSolver.hh
 *  @brief  LinearSolver class definition.
 *  @note   Copyright (C) 2013 Jeremy Roberts.
 */
//----------------------------------------------------------------------------//

#ifndef linear_algebra_LINEARSOLVER_HH_
#define linear_algebra_LINEARSOLVER_HH_

#include "Vector.hh"
#include "MatrixBase.hh"
#include "comm/Comm.hh"
#include "utilities/DBC.hh"
#include "utilities/SP.hh"
#include "petsc.h"
#include <vector>
#include <string>

namespace linear_algebra
{

class LinearSolver
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<LinearSolver>  SP_solver;
  typedef MatrixBase::SP_matrix               SP_matrix;
  typedef Vector::SP_vector                   SP_vector;

  //--------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param A  pointer to linear system matrix
   *  @param P  pointer to preconditioning matrix, possibly equal to A
   */
  LinearSolver(SP_matrix A, SP_matrix P);


  /**
   *  @brief Solve the linear system
   *  @param b  right hand side
   *  @param x  solution
   */
  void solve(SP_vector b, SP_vector x);

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// PETSc linear solver
  KSP d_solver;
  /// Matrix of linear system
  SP_matrix d_A;
  /// Preconditioner
  SP_matrix d_P;


};

} // end namespace linear_algebra

#endif /* linear_algebra_LINEARSOLVER_HH_ */

//----------------------------------------------------------------------------//
//              end of file LinearSolver.hh
//----------------------------------------------------------------------------//
