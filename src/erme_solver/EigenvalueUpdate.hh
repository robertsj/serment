//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  EigenvalueUpdate.hh
 *  @brief EigenvalueUpdate
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef erme_solver_EIGENVALUEUPDATE_HH_
#define erme_solver_EIGENVALUEUPDATE_HH_

#include "linear_algebra/Vector.hh"
#include "utilities/Definitions.hh"
#include "utilities/SP.hh"

namespace erme_solver
{

/**
 *  @class EigenvalueUpdate
 *  @brief Provides an updated estimate of the k-eigenvalue
 *
 *  The outer component of @ref GlobalSolverPicard involves updating the
 *  k-eigenvalue.  This is typically done via balance, but various techniques
 *  can be used to provide a better estimate, possibly via extrapolation.
 */
class EigenvalueUpdate
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<EigenvalueUpdate>  SP_update;
  typedef linear_algebra::Vector::SP_vector       SP_vector;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /// Constructor
  EigenvalueUpdate()
  {/* ... */}

  /// Virtual destructor
  virtual ~EigenvalueUpdate(){}

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /**
   *  @brief Computes and updated keff, possibly based on previous history
   *  @param keff   latest eigenvalue estimate
   *  @param J      latest boundary unknowns
   */
  virtual double compute(const double keff, SP_vector J)
  {
    /// Default does nothing
    return keff;
  }

};

} // end namespace erme_solver

#endif // erme_solver_EIGENVALUEUPDATE_HH_

//----------------------------------------------------------------------------//
//              end of file EigenvalueUpdate.hh
//----------------------------------------------------------------------------//
