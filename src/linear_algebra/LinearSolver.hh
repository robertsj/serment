//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LinearSolver.hh
 * \author robertsj
 * \date   Sep 5, 2012
 * \brief  LinearSolver class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef LINEARSOLVER_HH_
#define LINEARSOLVER_HH_

// Linear Algebra
#include "Vector.hh"
#include "MatrixBase.hh"

// Comm
#include "comm/Comm.hh"

// Detran Utilities
#include "DBC.hh"
#include "SP.hh"

// System
#include "petsc.h"
#include <vector>
#include <string>

namespace linear_algebra
{

class LinearSolver
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran::SP<LinearSolver>    SP_solver;
  typedef MatrixBase::SP_matrix       SP_matrix;
  typedef Vector::SP_vector           SP_vector;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   *  \param A  Pointer to linear system matrix
   *  \param P  Pointer to preconditioning matrix, possibly equal to A
   */
  LinearSolver(SP_matrix A, SP_matrix P);


  /*!
   *  \brief Solve the linear system
   *  \param b  right hand side
   *  \param x  solution
   */
  void solve(SP_vector b, SP_vector x);

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// PETSc linear solver
  KSP d_solver;

  /// Matrix of linear system
  SP_matrix d_A;

  /// Preconditioner
  SP_matrix d_P;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//


};

} // end namespace linear_algebra

#endif /* LINEARSOLVER_HH_ */

//---------------------------------------------------------------------------//
//              end of file LinearSolver.hh
//---------------------------------------------------------------------------//
