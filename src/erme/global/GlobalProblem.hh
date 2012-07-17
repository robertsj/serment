//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   GlobalProblem.hh
 * \author Jeremy Roberts
 * \date   11/23/2010
 * \brief  GlobalProblem class definition.
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 177                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-12-09 16:31:49 -0500 (Fri, 09 Dec 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef GLOBALPROBLEM_HH
#define GLOBALPROBLEM_HH

#include <iostream>
#include <cmath>
#include "LinAlg.hh"
#include "GlobalInput.hh"
#include "ResponseFunctionServer.hh"
#include "ResponseMatrix.hh"
#include "ResponseMatrixFull.hh"
#include "AbsorptionResponse.hh"
#include "FissionResponse.hh"
#include "LeakageResponse.hh"
#include "ConnectMatrix.hh"
#include "Connect2dCart.hh"

#include "utilities/SP.hh"

using namespace std;
//===========================================================================//
/*!
 * \class GlobalProblem
 * \brief Manages the global problem.
 *
 * to be completed
 *
 */
//===========================================================================//

class GlobalProblem
{

public:

  /// Typedefs
  //\{
  typedef typename util::SP<GlobalProblem>      SP_globalproblem;
  typedef typename GlobalInput::SP_globalinput  SP_globalinput;
  typedef typename ResponseMatrixFull::SP_R     SP_R;
  typedef typename Connect2dCart::SP_M          SP_M;
  //\}

  GlobalProblem(SP_globalinput input, ResponseFunctionServer *s)
    :  M(new Connect2dCart(input)),
       R(new ResponseMatrixFull(*input, s)),
       L(*input, s, M->getMindex(), M->getMsize()),
       F(*input, s),
       A(*input, s),
       J(input->degfree),
       d_run_already(false)
  {}

  ~GlobalProblem(){}

  /*!
   *  \brief Indicate whether the problem has been run before.
   *
   *  This is useful if the same problem is to be run repeatedly
   *  for timing studies, or if previous solutions should be used
   *  as the starting guess for subsequent cases.
   *
   */

  bool run_already() const { return d_run_already; }

  void set_run_already(bool run_already=true) { d_run_already = run_already; }

public:

  /// Connectivity matrix.  \todo Should not be 2d
  SP_M  M;

  /// Response matrix.
  SP_R  R;

  /// Leakage operator
  LeakageResponse L;

  /// Fission operator
  FissionResponse F;

  /// Absorption operator
  AbsorptionResponse A;

  /// Incident current vector
  SermentVector J;

private:

  /// Flag if this problem has been run before.
  bool d_run_already;

};

#endif // GLOBALPROBLEM_HH

//---------------------------------------------------------------------------//
//                 end of GlobalProblem.hh
//---------------------------------------------------------------------------//

