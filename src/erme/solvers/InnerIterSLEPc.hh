//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   InnerIterSLEPc.hh
 * \author Jeremy Roberts
 * \date   Nov 4, 2011
 * \brief  InnerIterSLEPC class definition.
 * \note   Copyright (C) 2011 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//
// $Rev::                                               $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date::                                              $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef INNERITERSLEPC_HH_
#define INNERITERSLEPC_HH_

#include "serment_config.h"

#include "slepceps.h"

#include "InnerIterBase.hh"

//===========================================================================//
/*!
 * \class InnerIterSLEPc
 * \brief Performs inner iterations using SLEPc within outer iteration.
 *
 * This is an interface class for use of SLEPc eigensolvers for inner
 * iterations.  SLEPc offers several algorithms which can be further
 * supplemented by additional libraries (e.g. ARPACK).
 *
 */
//===========================================================================//

class InnerIterSLEPc : public InnerIterBase
{


public:

  /// Typedefs
  //\{
  typedef typename util::SP<InnerIterSLEPc>     SP_inneriter;
  //\}

  /*!
   *  \brief Constructor
   *
   */
  InnerIterSLEPc(SP_M M, SP_R R);

  /*!
   *  \brief Destructor.
   *
   *  Frees the SLEPc items.
   *
   */
  ~InnerIterSLEPc();


  /*!
   *  \brief Perform inner iterations.
   *
   *  \param    max_iterations      Maximum number of iterations, nominally
   *                                meaning applications of MR on a vector.
   *  \param    tolerance           Generic tolerance on convergence.
   *  \param    X_0                 Initial guess.
   */
  scalar solve(int max_iters, double tol, SP_vector J_in, SP_vector J);

  // Declare the matrix shell wrapper a friend.
  friend inline PetscErrorCode apply_MR(Mat A, Vec V_in, Vec V_out);

private:

  /// SLEPc eigensolver context
  EPS       d_eps;

  /// SLEPc spectral transformation context
  ST        d_st;

};

//---------------------------------------------------------------------------//
/*!
 * \brief A matrix-vector multiplication wrapper function for MR.
 *
 * This is needed because PETSc works with function pointers, which does not
 * include pointers to member functions.
 *
 * \param   A       PETSc shell matrix
 * \param   X_in    Incoming PETSc vector
 * \param   X_out   Outgoing PETSc vector
 *
 */
inline PetscErrorCode apply_MR(Mat A, Vec X_in, Vec X_out);

#endif /* INNERITERSLEPC_HH_ */

//---------------------------------------------------------------------------//
//              end of InnerIterSLEPc.hh
//---------------------------------------------------------------------------//
