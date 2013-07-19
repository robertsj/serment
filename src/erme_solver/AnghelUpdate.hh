//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  AnghelUpdate.hh
 *  @brief AnghelUpdate
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//


#ifndef erme_solver_ANGHELUPDATE_HH_
#define erme_solver_ANGHELUPDATE_HH_

#include "EigenvalueUpdate.hh"
#include "utilities/DBC.hh"
#include <cmath>
#include <cstdio>

namespace erme_solver
{

/**
 *  @class AnghelUpdate
 *  @brief Update based on interpolation using an exponential form
 *
 *  Anghel and Gheorghiu state that current responses follow an exponential
 *  dependence on @f$ k @f$ and assume that so does @f$ \lambda @f$, i.e.
 *  @f[
 *    \lambda(k) \approx a e^{b/k} \, .
 *  @f]
 *  Using an initial guess @f$ k_0 @f$ (with corresponding  @f$ \lambda_0 @f$
 *  and updating for @f$ k_1 @f$ and @f$ \lambda_1 @f$ via balance, one can
 *  define the coefficients @f$ a @f$ and @f$ b @f$.  Then, one can solve for
 *  @f$ k_2 @f$ such that @f$ \lambda_2 = 1 @f$, since, when the responses are
 *  conservative, the critical @f$ k @f$ leads to @f$ \lambda = 1 @f$.
 *  Of course, this scheme is problematic if the approximations or
 *  convergence criteria fail to yield @f$ \lambda = 1 @f$.
 *
 *  Actually, similar schemes have been used earlier.  For instance, in his
 *  thesis, Lindahl uses an interpolation based on the assumed relation
 *  @f[
 *    \frac{1}{k} \propto \frac{1}{\lambda}
 *  @f]
 *  coupled with Regula Falsi (linear extrapolation). Beyond the third
 *  estimate, a parabolic extrapolation is used.
 *  Forget (unfamiliar at the time with Lindah's work) has a recent summary
 *  that uses
 *  @f[
 *    k \propto \frac{1}{\lambda} \, .
 *  @f]
 *
 *  Here, we implement a two-term interpolation based on the exponential form
 *  as well as the forms @f$ k \propto \lambda @f$,
 *  @f$ k \propto \lambda^{-1} @f$, and @f$ k^{-1} \propto \lambda^{-1} @f$.
 *  Each is based on the assumption that @f$ \lambda = 1 @f$ tends to
 *  exactly unity.
 *
 *  See Anghel and Gheorghiu, "An Iteration Strategy for the Response
 *    Matrix Method," Ann. Nucl. Energy, <em>14</em, 5, pp 219-226 (1986)
 */
class AnghelUpdate: public EigenvalueUpdate
{

public:

  //--------------------------------------------------------------------------//
  // ENUMERATIONS
  //--------------------------------------------------------------------------//

  enum INTERPOLATION_SCHEMES
  {
    EXP, LINLIN, INVLIN, INVINV, END_INTERPOLATION_SCHEMES
  };
  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::vec_dbl  vec_dbl;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /// Constructor
  AnghelUpdate(const int scheme = EXP, bool keep_new = true)
    : d_scheme(scheme)
    , d_counter(0)
    , d_keep_new(keep_new)
  {
    Require(scheme < END_INTERPOLATION_SCHEMES);
    // Assumed initial guess
    d_keff[0] = 1.0;
  }

  /// Virtual destructor
  virtual ~AnghelUpdate(){}

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// Computes an updated keff, possibly based on previous history
  double compute(const double keff, const double lambda, SP_vector J)
  {
    double new_keff = keff;
    if (d_counter == 0)
    {
      d_keff[1]   = keff;
      d_lambda[0] = lambda;
    }
    else
    {
      if (d_counter > 1) d_lambda[0] = d_lambda[1];
      d_lambda[1] = lambda;
      new_keff = update(keff, d_keff[0], d_keff[1], d_lambda[0], d_lambda[1]);
      d_keff[2] = d_keep_new ? new_keff : keff;
      d_keff[0] = d_keff[1];
      d_keff[1] = d_keff[2];
    }
    ++d_counter;
    return new_keff;
  }

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Interpolation scheme
  int d_scheme;
  /// Iteration counter
  int d_counter;
  /// Values of @f$ k @f$ and  @f$ \lambda @f$ for past two iterations
  //@{
  double d_keff[3];
  double d_lambda[2];
  //@}
  /// Keep updated keff rather than the inital unadjusted keff as last keff
  bool d_keep_new;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  double update(const double korig,
                const double k0,
                const double k1,
                const double l0,
                const double l1)
  {
    Require(k0 > 0);
    Require(k1 > 0);
    Require(l0 > 0);
    Require(l1 > 0);
//    Require(k0 != k1);
//    Require(l0 != l1);

    double keff = korig, a = 0.0, b = 0.0;
    if (d_scheme == EXP)
    {
      a = l0 * std::pow(l0 / l1, k1 / (k0 - k1));
      b = std::log(l0 / l1) * k1 * k0 / (k1 - k0);
      double bb = -1;
      if (std::abs(l0 - l1) < 1e-9)
      {
        // small argument approximation for log
        double X = (l0 - l1) / l1;
        X = X - 0.5 * X * X + (1.0/3.0) * X * X * X;
        bb = X * k1 * k0 / (k1 - k0);
      }
      // Anghel suggests these criteria indicate a good fit.  Otherwise, it's
      // subject to blow up.  In that case, stick with the k from balance.
      if (a < 1.0 && b > 0.0)
      {
        //std::printf(" %20.16e  %20.16e   %20.16e \n", b, bb, a);
        if (bb > 0) b = bb;
        keff = -b / std::log(a);
      }
      else
      {
        std::cout << " ouch " << a << std::endl;
      }
    }
    else if (d_scheme == LINLIN)
    {
      if (std::abs(k0 - k1) < 1e-10 || std::abs(l0 - l1) < 1e-10)
      {
        // avoid roundoff
        keff = (1.0-l0)*(k1-k0)/(l1-l0)+k0;
      }
      else
      {
        a = (l0 - l1) / (k0 - k1);
        b = (k0 * l1 - k1 * l0) / (k0 - k1);
        keff = (1.0 - b) / a;
      }
    }
    else if (d_scheme == INVLIN)
    {
      if (std::abs(k0 - k1) < 1e-9 || std::abs(l0 - l1) < 1e-9)
      {
        // avoid roundoff
        a = (1.0 - l0) * (k1 - k0);
        keff = a / (l1 - l0) + k0;
      }
      else
      {
        a = -k0 * k1 * (l0 - l1) / (k0 - k1);
        b = (k0 * l0 - k1 * l1) / (k0 - k1);
        keff = a / (1.0 - b);
      }
    }
    else if (d_scheme == INVINV)
    {
      if (std::abs(k0-k1) < 1e-8 || std::abs(l0-l1) < 1e-8)
      {
        // avoid roundoff
        a = (1.0 - 1.0 * l0) * (k1 - k0);
        keff = a * (1.0 + l0 / (l1 - l0)) + k0;
      }
      else
      {
        a = k0 * k1 * (l0 - l1) / (l0 * l1 * (k0 - k1));
        b = (k0 * l1 - k1 * l0) / (l0 * l1 * (k0 - k1));
        keff = -a/(b-1);
      }
    }
    return keff;
  }

};

} // end namespace erme_solver

#endif /* erme_solver_ANGHELUPDATE_HH_ */

//----------------------------------------------------------------------------//
//              end of file AnghelUpdate.hh
//----------------------------------------------------------------------------//
