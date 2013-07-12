//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  NonlinearResidualBase.hh
 *  @brief NonlinearResidualBase
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//


#ifndef linear_algebra_NONLINEARRESIDUALBASE_HH_
#define linear_algebra_NONLINEARRESIDUALBASE_HH_

#include "Vector.hh"
#include "utilities/Definitions.hh"
#include "utilities/SP.hh"

namespace linear_algebra
{

/// Base class for nonlinear residuals
class NonlinearResidualBase
{

public :

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<NonlinearResidualBase>    SP_residual;
  typedef detran_utilities::size_t                       size_t;

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// Virtual destructor
  virtual ~NonlinearResidualBase(){}

  /// Evaluate the residual vector given an approximate solution x
  virtual void evaluate(Vector &x, Vector &f) = 0;

};

} // end namespace linear_algebra


//----------------------------------------------------------------------------//
//                 end of NonlinearResidualBase.hh
//----------------------------------------------------------------------------//


#endif /* linear_algebra_NONLINEARRESIDUALBASE_HH_ */
