//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   SteffensenUpdate.hh
 *  @brief  SteffensenUpdate
 *  @author Jeremy Roberts
 *  @date   Oct 4, 2012
 */
//---------------------------------------------------------------------------//

#ifndef erme_solver_STEFFENSENUPDATE_HH_
#define erme_solver_STEFFENSENUPDATE_HH_

#include "EigenvalueUpdate.hh"

namespace erme_solver
{

class SteffensenUpdate: public EigenvalueUpdate
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::vec_dbl  vec_dbl;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /// Constructor
  SteffensenUpdate()
    : d_minimum(4)
  { /* ... */ }

  /// Virtual destructor
  virtual ~SteffensenUpdate(){}

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /**
   *  @brief Computes and updated keff, possibly based on previous history
   *  @param keff   latest eigenvalue estimate
   *  @param J      latest boundary unknowns
   */
  double compute(const double keff, SP_vector J)
  {
    // Save the latest keff
    d_keffs.push_back(keff);
    // Check size
    int n = d_keffs.size();
    if (d_keffs.size() < d_minimum) return keff;
    // If sequence is long enough, extrapolate
    double k0 = d_keffs[n-3];
    double k1 = d_keffs[n-2];
    double k2 = d_keffs[n-1];
    return k0 - (k1 - k0)*(k1 - k0) / (k2 - 2.0*k1 + k0);
  }

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Values of keff for each iteration
  vec_dbl d_keffs;
  /// Minimum number of computed keffs before extrapolating
  int d_minimum;

};

} // end namespace erme_solver

#endif // erme_solver_STEFFENSENUPDATE_HH_

//---------------------------------------------------------------------------//
//              end of file SteffensenUpdate.hh
//---------------------------------------------------------------------------//
