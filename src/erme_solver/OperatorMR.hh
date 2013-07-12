//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   OperatorMR.hh
 *  @brief  OperatorMR class definition
 *  @note   Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

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

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<OperatorMR>          SP_MR;
  typedef erme::ResponseMatrix::SP_responsematrix   SP_R;
  typedef erme::Connect::SP_connect                 SP_M;
  typedef linear_algebra::Vector                    Vector;

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

  void multiply(Vector &v_in, Vector &v_out)
  {
     d_R->multiply(v_in, d_V);
     d_M->multiply(d_V, v_out);
  }

  void multiply_transpose(Vector &v_in, Vector &v_out)
  {
    d_R->multiply_transpose(v_in, d_V);
    d_M->multiply_transpose(d_V, v_out);
  }

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Response matrix
  SP_R d_R;
  /// Connectivity matrix
  SP_M d_M;
  /// Temporary vector for intermediate R*X
  linear_algebra::Vector d_V;

};


} // end namespace erme_solver

#endif // erme_solver_OPERATORMR_HH_

//----------------------------------------------------------------------------//
//              end of file OperatorMR.hh
//----------------------------------------------------------------------------//
