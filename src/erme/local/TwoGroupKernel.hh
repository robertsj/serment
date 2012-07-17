//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   TwoGroupKernel.hh
 * \author Jeremy Roberts
 * \date   10/14/2011
 * \brief  TwoGroupKernel class definition.
 * \note   Copyright (C) 2011 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 140                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-09-14 12:53:40 -0400 (Wed, 14 Sep 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#include "Diff2dInput.hh"
#include "Diff2dProblem.hh"
#include "Diff2dSolver.hh"
#include "ResponseFunction.hh"
#include "utilities/SP.hh"


namespace erme
{

//===========================================================================//
/*!
 * \class TwoGroupKernel
 * \brief Simple two-dimensional, two-group response function generator.
 *
 * This class produces response functions for homogeneous cells based on a
 * semi-analytical model augmented with fitting parameters.  For the tests
 * used to guide fitting, the response functions are all within about 10%
 * relative error, and in absolute terms, the residuals were usually below
 * 0.001.
 *
 */
//===========================================================================//

class TwoGroupKernel
{

public:

  /// \name Typedefs
  //\{
  typedef Diff2dInput::SP_input             SP_input;
  typedef Diff2dProblem::SP_problem         SP_problem;
  typedef Diff2dSolver::SP_solver           SP_solver;
  typedef util::SP<ResponseFunction>        SP_responsefunction;
  typedef std::vector<double>               Vec_Dbl;
  typedef std::vector<Vec_Dbl>              Mat_Dbl;
  //\}

  /*!
   *  \brief Constructor
   *
   */
  TwoGroupKernel();

  /*!
   *  \brief Destructor
   *
   */
  ~TwoGroupKernel();

  /*!
   *  \brief Get the responses.
   *
   *  Responses are returned as a function of the following:
   *  + Delta  -- Assembly dimension
   *  + k      -- k-effective
   *  + D1     -- Group 1 diffusion coefficient
   *  + D2     -- Group 2 diffusion coefficient
   *  + R1     -- Group 1 removal cross-section
   *  + A2     -- Group 2 absorption cross-section
   *  + F1     -- Group 1 fission cross-section times nu
   *  + F2     -- Group 2 fission cross-section times nu
   *  + S12    -- Group 1 to 2 scattering cross-section
   *  These data are held sequentially in the data vector.
   *
   *  \param  responses     Response function smart pointer to be filled.
   *  \param  data          Two group data.
   */
  void response(SP_responsefunction responses, Vec_Dbl &data);


private:

  /// \name ???
  //\{

  //\}

  /// \name Fit Data
  //\{

  /// Fit parameters
  double beta1[40];
  double beta2[40];

  //\}
};




} // end namespace erme

//---------------------------------------------------------------------------//
//                 end of TwoGroupKernel.hh
//---------------------------------------------------------------------------//

