//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  NonlinearSolver.hh
 *  @brief NonlinearSolver class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef linear_algebra_NONLINEARSOLVER_HH_
#define linear_algebra_NONLINEARSOLVER_HH_

#include "JacobianBase.hh"
#include "NonlinearResidualBase.hh"
#include "Vector.hh"
#include "comm/Comm.hh"
#include "utilities/InputDB.hh"
#include "utilities/SP.hh"
#include "petsc.h"
#include <vector>
#include <string>

namespace linear_algebra
{

/**
 *  @class NonlinearSolver
 *  @brief Lightweight wrapper arond SNES
 *
 *  The required SNES functionality is encapsulated.  Relevant database
 *  entries include
 *    -- newton_type
 *       0 : Standard Newton method based on a Jacobian matrix given
 *           by the user.  This can include Jacobian shell matrices in which
 *           only the action is given.
 *       1 : Jacobian-Free Newton-Krylov via finite-differencing the
 *           residual function.
 *       2 : Newton via a full Jacobian constructed from finite-differencing
 *           the residual.  Useful primarily for testing.
 */
class NonlinearSolver
{

public:

  //--------------------------------------------------------------------------//
  // ENUMERATIONS
  //--------------------------------------------------------------------------//

  /**
   *  NEWTON -- regular newton based on whatever Jacobian matrix is given
   *  JFNK   -- approximate Jacobian action via finite difference of f(x)
   *  FDJAC  -- full Jacobian via finite difference (for debug only)
   */
  enum newton_types
  {
    NEWTON, JFNK, FDJAC, END_NEWTON_TYPES
  };

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<NonlinearSolver>       SP_solver;
  typedef JacobianBase::SP_jacobian                   SP_jacobian;
  typedef MatrixBase::SP_matrix                       SP_matrix;
  typedef Vector::SP_vector                           SP_vector;
  typedef NonlinearResidualBase::SP_residual          SP_residual;
  typedef detran_utilities::InputDB::SP_input         SP_db;

  //--------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //--------------------------------------------------------------------------//

  /// Constructor
  NonlinearSolver();

  /**
   *  @brief Setup the solver
   *  @param   db  parameter database
   *  @param    J  jacobian matrix
   *  @param    P  preconditioner matrix (may be J)
   */
  void setup(SP_db db, SP_residual f, SP_jacobian J, SP_jacobian P);

  /// Solver the nonlinear system for an initial guess x
  void solve(SP_vector    x,
             const double atol = 1e-8,
             const double rtol = 1e-8,
             const int    maxit = 20);

  //@{
  ///  Getters
  int number_iterations() const;
  int number_linear_iterations() const;
  SP_residual residual();
  SP_jacobian jacobian();
  SP_jacobian preconditioner();
  //@}

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Matrix of linear system
  SP_jacobian d_J;
  /// Preconditioner
  SP_jacobian d_P;
  /// Residual function
  SP_residual d_residual;
  /// The residual vector
  SP_vector d_r;
  /// PETSc linear solver
  SNES d_solver;
  /// Number of iterations
  int d_number_iterations;
  /// Number of linear solver iterations
  int d_number_linear_iterations;
  /// Temporary matrix for JFNK Jacobian
  SP_matrix d_jfnk_jacobian;

};

//----------------------------------------------------------------------------//
inline PetscErrorCode
residual_wrap(SNES snes, Vec x, Vec f, void *context)
{
  //std::cout << " residual_wrap WRAP" << std::endl;

  PetscErrorCode ierr = 0;
  NonlinearSolver* foo = (NonlinearSolver*) context;
  Vector X(x), F(f);
  //std::cout << " X =" << std::endl;
  //X.display();
  //std::cout << " end X" << std::endl;
  foo->residual()->evaluate(&X, &F);
  return  ierr;
};

//----------------------------------------------------------------------------//
inline PetscErrorCode
jacobian_wrap(SNES snes, Vec x, Mat *J, Mat *P, MatStructure *flag, void *ctx)
{
  //std::cout << " jacobian_wrap WRAP" << std::endl;

  PetscErrorCode ierr = 0;
  NonlinearSolver* foo = (NonlinearSolver*) ctx;
  Vector::SP_vector X(new Vector(x));
  // update the jacobian matrix
  foo->jacobian()->update(X);
  // if we have two distinct matrices, update the second
  if (J != P) foo->preconditioner()->update(X);
  return ierr;
};

//----------------------------------------------------------------------------//
/**
 *  This is the function PETSc calls to update the Jacobian and its
 *  preconditioner when using JFNK.  That is, the action of the Jacobian
 *  is approximated by a finite difference (and so the matrix passed
 *  needs to update beyond assembly).  However, the preconditioner can still
 *  be updated, as done here.
 */
inline PetscErrorCode
jfnk_wrap(SNES snes, Vec x, Mat *J, Mat *P, MatStructure *flag, void *ctx)
{
  PetscErrorCode ierr = 0;
  //std::cout << " JFNK WRAP" << std::endl;
  // Assemble the Jacobian matrix (this is all MatMFFDComputeJacobian does)
  ierr = MatAssemblyBegin(*J, MAT_FINAL_ASSEMBLY);
  ierr = MatAssemblyEnd(*J, MAT_FINAL_ASSEMBLY);

  // Update the preconditioner P, if distinct from J
  if (J != P)
  {
    NonlinearSolver* foo = (NonlinearSolver*) ctx;
    Vector::SP_vector X(new Vector(x));
    foo->preconditioner()->update(X);
  }

  return ierr;
}

} // end namespace linear_algebra

#endif /* linear_algebra_NONLINEARSOLVER_HH_ */

//----------------------------------------------------------------------------//
//                 end of NonlinearSolver.hh
//----------------------------------------------------------------------------//
