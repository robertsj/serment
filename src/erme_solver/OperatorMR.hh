//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   OperatorMR.hh
 *  @brief  OperatorMR class definition
 *  @author Jeremy Roberts
 *  @date   Oct 4, 2012
 */
//---------------------------------------------------------------------------//

#ifndef erme_solver_OPERATORMR_HH_
#define erme_solver_OPERATORMR_HH_

#include "erme/Connect.hh"
#include "erme/ResponseMatrix.hh"
#include "linear_algebra/MatrixShell.hh"

namespace erme_solver
{

/**
 *  @class OperatorMR
 *  @brief Performs the action of M*R
 */
class OperatorMR: public linear_algebra::MatrixShell
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<OperatorMR>          SP_MR;
  typedef erme::ResponseMatrix::SP_responsematrix   SP_R;
  typedef erme::Connect::SP_connect                 SP_M;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param R  Response matrix
   *  @param M  Connectivity matrix
   */
  OperatorMR(SP_R R, SP_M M)
    : MatrixShell(R->number_local_rows(), R->number_local_columns(), this)
    , d_R(R)
    , d_M(M)
    , d_V(R->number_local_columns(), 0.0)
  {
    Require(d_R);
    Require(d_M);
  }

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Response matrix
  SP_R d_R;
  /// Connectivity matrix
  SP_M d_M;
  /// Temporary vector for intermediate R*X
  linear_algebra::Vector d_V;

  //---------------------------------------------------------------------------//
  // SHELL MATRIX OPERATIONS
  //---------------------------------------------------------------------------//

  /**
   *  @brief Matrix-vector multiplication
   *  @param x  Input vector
   *  @param y  Output vector
   */
  PetscErrorCode shell_multiply(Vec x, Vec y)
  {
    PetscErrorCode ierr;
    // Compute R*x -> v
    ierr = MatMult(d_R->A(), x, d_V.V());
    Ensure(!ierr);
    // Compute M*v -> y
    ierr = MatMult(d_M->A(), d_V.V(), y);
    Ensure(!ierr);
    return ierr;
  }

  /**
   *  @brief Matrix-vector multiplication using matrix transpose.
   *  @param x  Input vector
   *  @param y  Output vector
   */
  PetscErrorCode shell_multiply_transpose(Vec x, Vec y)
  {
    // (Note that MR)' = R'M'
    PetscErrorCode ierr;
    // Compute M'*x -> v
    ierr = MatMultTranspose(d_M->A(), x, d_V.V());
    Ensure(!ierr);
    // Compute R'*v -> y
    ierr = MatMultTranspose(d_R->A(), d_V.V(), y);
    Ensure(!ierr);
    return ierr;
  }

};


} // end namespace erme_solver

#endif // erme_solver_OPERATORMR_HH_

//---------------------------------------------------------------------------//
//              end of file OperatorMR.hh
//---------------------------------------------------------------------------//
