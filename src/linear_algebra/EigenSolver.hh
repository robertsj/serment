//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   EigenSolver.hh
 * \author robertsj
 * \date   Sep 5, 2012
 * \brief  EigenSolver class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef EIGENSOLVER_HH_
#define EIGENSOLVER_HH_

// Linear Algebra
#include "Vector.hh"
#include "MatrixBase.hh"

// Comm
#include "comm/Comm.hh"

// Detran Utilities
#include "DBC.hh"
#include "SP.hh"

// System
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

  typedef detran::SP<EigenSolver>    SP_solver;
  typedef MatrixBase::SP_matrix       SP_matrix;
  typedef Vector::SP_vector           SP_vector;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   *  \param A  Pointer to left hand operator
   *  \param P  Pointer to right hand operator (possibly null)
   */
  EigenSolver(SP_matrix A, SP_matrix B = SP_matrix(0));


  /*!
   *  \brief Solve the eigenvalue problem
   *  \param x  Eigenvector
   */
  void solve(SP_vector x);

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

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//


};

#endif /* EIGENSOLVER_HH_ */

//---------------------------------------------------------------------------//
//              end of file EigenSolver.hh
//---------------------------------------------------------------------------//
