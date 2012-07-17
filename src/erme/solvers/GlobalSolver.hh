//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   GlobalSolver.hh
 * \author Jeremy Roberts
 * \date   11/23/2010
 * \brief  GlobalSolver class definition.
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 167                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-11-07 23:01:04 -0500 (Mon, 07 Nov 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef GLOBALSOLVER_HH
#define GLOBALSOLVER_HH

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
#include "GlobalProblem.hh"

#include "utilities/SP.hh"

using namespace std;
//===========================================================================//
/*!
 * \class GlobalSolver
 * \brief A base class for all global solvers.
 *
 * The eigenvalue
 *
 */
//===========================================================================//

class GlobalSolver
{

public:

  /// Typedefs
  //\{
  typedef typename util::SP<GlobalSolver>           SP_globalsolver;
  typedef typename GlobalProblem::SP_globalproblem  SP_globalproblem;
  typedef typename GlobalInput::SP_globalinput      SP_globalinput;
  typedef typename SermentVector::SP_vector         SP_vector;
  //\}

  /*!
   *  \brief Constructor
   *
   */
  GlobalSolver(SP_globalproblem problem, SP_globalinput input);

  /*!
   *  \brief Pure virtual destructor.
   *
   */
  ~GlobalSolver();

  /*!
   *  \brief Solve the eigenvalue response matrix equations.
   *
   *
   */
  virtual void solve() = 0;

protected:

  /// \name Data
  //\{

  /// Global problem specification.
  SP_globalproblem      d_problem;

  /// Global input data.
  SP_globalinput        d_input;

  ///
  vector<scalar>        keffHistory;

  ///
  vector<scalar>        keffErrHistory;

  ///
  vector<scalar>        curErrHistory;

  ///
  vector<scalar>        residErrHistory;

  ///
  scalar                wallTime;

  //\}

  /// \name Implementation
  //\{

  /*!
   *  \brief Compute coarse mesh fission rates.
   *
   *  Computes \f$ \nu\Sigma_f\phi \f$ in each mesh.
   */
  void fissionRates(SP_vector J);

  //\}

};

#endif // GLOBALSOLVER_HH
//---------------------------------------------------------------------------//
//                 end of GlobalSolver.hh
//---------------------------------------------------------------------------//

