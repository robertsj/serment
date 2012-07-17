//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   InnerPowerIter.hh
 * \author Jeremy Roberts
 * \date   Sep 24, 2011
 * \brief  InnerIterPower class definition.
 * \note   Copyright (C) 2011 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//
// $Rev::                                               $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date::                                              $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef INNERITERPOWER_HH_
#define INNERITERPOWER_HH_

#include "InnerIterBase.hh"

//===========================================================================//
/*!
 * \class InnerIterPower
 * \brief Performs inner iterations via power method within outer iteration.
 *
 * This is our own implementation of the power method for inners within
 * power-like outers.  The SLEPc interface also provides the power method,
 * and their implementation is to be preferred.
 *
 * The power (iteration) method is a standard procedure for finding the
 * largest eigenvalue of an operator.
 *
 *
 */
//===========================================================================//

class InnerIterPower : public InnerIterBase
{


public:

  /// Typedefs
  //\{
  typedef typename util::SP<InnerIterPower>     SP_inneriter;
  //\}

  /*!
   *  \brief Constructor
   *
   */
  InnerIterPower(SP_M M, SP_R R);

  /*!
   *  \brief Destructor
   *
   */
  ~InnerIterPower(){}


  /*!
   *  \brief Perform inner iterations.
   *
   *  \param    max_iterations      Maximum number of iterations, nominally
   *                                meaning applications of MR on a vector.
   *  \param    tolerance           Generic tolerance on convergence.
   *  \param    X_0                 Initial guess.
   */
  scalar solve(int max_iters, double tol, SP_vector J_in, SP_vector J);

};

#endif /* INNERITERPOWER_HH_ */

//---------------------------------------------------------------------------//
//              end of InnerIterPower.hh
//---------------------------------------------------------------------------//
