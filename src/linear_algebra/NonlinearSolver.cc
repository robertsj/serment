//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  NonlinearSolver.cc
 *  @brief NonlinearSolver
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "NonlinearSolver.hh"
#include "Matrix.hh"
#include "utilities/DBC.hh"

namespace linear_algebra
{

//----------------------------------------------------------------------------//
NonlinearSolver::NonlinearSolver()
  : d_number_iterations(0)
  , d_number_linear_iterations(0)
{
  SNESCreate(PETSC_COMM_WORLD, &d_solver);
}

//----------------------------------------------------------------------------//
void NonlinearSolver::setup(SP_db       db,
                            SP_residual f,
                            SP_jacobian J,
                            SP_jacobian P)
{
  Require(db);
  Require(J);
  Require(P);
  d_J = J;
  d_P = P;
  d_residual = f;
  d_r = new Vector(d_J->local_size(), 0.0);

  // set the residual function
  SNESSetFunction(d_solver, d_r->V(), residual_wrap, this);

  // newton type
  int newton_type = NEWTON;
  if (db->check("newton_type"))
    newton_type = db->get<int>("newton_type");
  Assert(newton_type < END_NEWTON_TYPES);
  if (newton_type == NEWTON)
  {

    SNESSetJacobian(d_solver,
                    d_J->matrix()->A(), d_P->matrix()->A(),
                    jacobian_wrap, this);
  }
  else if (newton_type == JFNK)
  {
    typedef PetscErrorCode (*compute_f)(void*, Vec, Vec);
    d_jfnk_jacobian = new Matrix();
    Mat A = d_jfnk_jacobian->A();
    // this creates a "matrix free" matrix for finite-difference action
    MatCreateSNESMF(d_solver, &A);
    // this sets the function used in finite difference to call whatever
    // the residual is set to be (REDUNTANT; by default uses residual)
    //MatMFFDSetFunction(A, (compute_f) SNESComputeFunction, d_solver);
    // indicates that
    SNESSetJacobian(d_solver, A, d_P->matrix()->A(), jfnk_wrap, this);
  }
  else if (newton_type == FDJAC)
  {
    // this computes an explicit jacobian using finite differences
    SNESSetJacobian(d_solver, d_J->matrix()->A(), d_P->matrix()->A(),
                    SNESDefaultComputeJacobian, NULL);
  }

  SNESSetFromOptions(d_solver);
}

//----------------------------------------------------------------------------//
void NonlinearSolver::solve(SP_vector x)
{
  Require(x);
  PetscErrorCode ierr = 0;
  ierr = SNESSolve(d_solver, PETSC_NULL, x->V());
  ierr = SNESGetIterationNumber(d_solver, &d_number_iterations);
  ierr = SNESGetLinearSolveIterations(d_solver, &d_number_linear_iterations);
  Ensure(!ierr);
}

//----------------------------------------------------------------------------//
int NonlinearSolver::number_iterations() const
{
  return d_number_iterations;
}

//----------------------------------------------------------------------------//
int NonlinearSolver::number_linear_iterations() const
{
  return d_number_linear_iterations;
}

//----------------------------------------------------------------------------//
NonlinearSolver::SP_residual NonlinearSolver::residual()
{
  return d_residual;
}

//----------------------------------------------------------------------------//
NonlinearSolver::SP_jacobian NonlinearSolver::jacobian()
{
  Require(d_J);
  return d_J;
}

//----------------------------------------------------------------------------//

NonlinearSolver::SP_jacobian NonlinearSolver::preconditioner()
{
  Require(d_P);
  return d_P;
}

} // end namespace linear_algebra

//----------------------------------------------------------------------------//
//                 end of NonlinearSolver.cc
//----------------------------------------------------------------------------//
