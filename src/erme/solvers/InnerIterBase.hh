//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   InnerIterBase.hh
 * \author Jeremy Roberts
 * \date   Nov 4, 2011
 * \brief  InnerIterBase class definition.
 * \note   Copyright (C) 2011 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//
// $Rev::                                               $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date::                                              $:Date of last commit
//---------------------------------------------------------------------------//


#ifndef INNERITERBASE_HH_
#define INNERITERBASE_HH_

#include "ResponseMatrixFull.hh"
#include "Connect2dCart.hh"

#include "linalg/LinAlg.hh"
#include "utilities/SP.hh"

//===========================================================================//
/*!
 * \class InnerIterBase
 * \brief Base class for inner iterations within an outer power-like iteration.
 */
//===========================================================================//
class InnerIterBase
{

public:

  /// Typedefs
  //\{
  typedef typename util::SP<InnerIterBase>      SP_inneriter;
  typedef typename ResponseMatrixFull::SP_R     SP_R;
  typedef typename Connect2dCart::SP_M          SP_M;
  typedef typename SermentVector::SP_vector     SP_vector;
  //\}

  /*!
   *  \brief Constructor
   *
   */
  InnerIterBase(SP_M M, SP_R R) :
    d_M(M), d_R(R), d_J_tmp(new SermentVector(R->m))
  {
  }

  virtual ~InnerIterBase()
  {
    //d_J_tmp->releaseMe();
  }


  /*!
   *  \brief Perform inner iterations.
   *
   *  \param    max_iterations      Maximum number of iterations, nominally
   *                                meaning applications of MR on a vector.
   *  \param    tolerance           Generic tolerance on convergence.
   *  \param    J_0                 Initial guess.
   *  \param    J                   Solved eigenvector.
   *  \return                       Dominant eigenvalue.
   */
  virtual scalar solve(int max_iters, double tol, SP_vector J_in, SP_vector J) = 0;

protected:

  /// \name Private Data
  //\{

  /// Maximum iterations
  int d_max_iters;

  /// Tolerance
  double d_tol;

  /// Connectivity Matrix
  SP_M d_M;

  /// Response Matrix
  SP_R d_R;

  /// Temporary storage for R*X_in = X_out
  SP_vector d_J_tmp;

  /// Shell matrix for use with mat-vec wrapper.
  Mat d_A;

  //\}

  /// \name Implementation
  //\{

  /*!
   *  \brief Return smart pointer to R matrix.
   *
   */
  SP_R R()
  {
    return d_R;
  }

  /*!
   *  \brief Return smart pointer to M matrix.
   *
   */
  SP_M M()
  {
    return d_M;
  }

  /*!
   *  \brief Mutable access to temporary vector.
   *
   */
  SP_vector J_tmp()
  {
    return d_J_tmp;
  }

  //\}


};

#endif /* INNERITERBASE_HH_ */

//---------------------------------------------------------------------------//
//              end of InnerIterBase.hh
//---------------------------------------------------------------------------//
