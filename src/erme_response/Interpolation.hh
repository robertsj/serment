//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Interpolation.hh
 *  @brief  Interpolation utilities
 *  @author Jeremy Roberts
 *  @date   Oct 9, 2012
 */
//---------------------------------------------------------------------------//

#ifndef erme_response_INTERPOLATION_HH_
#define erme_response_INTERPOLATION_HH_

#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"

namespace erme_response
{

/**
 *  @brief Interpolate a function evaluated at two points using
 *         linear interpolation
 */
inline double
interpolate_linear(double x, double x0, double x1, double r0, double r1)
{
//  double r = (r1 - r0) * (x - x0) / (x1 - x0) + r0;
//  std::cout << " x0 x1 r0 r1 = " << x0 << " " << x1 << " " << r0 << " " << r1 << " " << x << " " << x << " " << r << std::endl;
  return  (r1 - r0) * (x - x0) / (x1 - x0) + r0;
}

/**
 *  @brief Interpolate a function evaluated at three points using
 *         quadratic interpolation
 *
 *  This was generated using Maple's code generator with
 *  optimization.
 */
inline double
interpolate_quadratic(double x, double x0, double x1, double x2, double r0, double r1, double r2)
{
  double t10;
  double t11;
  double t13;
  double t14;
  double t15;
  double t16;
  double t27;
  double t3;
  double t33;
  double t43;
  double t6;
  double t8;
  t3 = x1 * r0;
  t6 = x0 * r1;
  t8 = x2 * x2;
  t10 = x0 * x0;
  t11 = x1 * t10;
  t13 = x2 * t10;
  t14 = x1 * x1;
  t15 = x2 * t14;
  t16 = x0 * t14;
  t27 = 0.1e1 / (x0 - x1);
  t33 = 0.1e1 / (-x0 * x2 + x0 * x1 + t8 - x1 * x2);
  t43 = x * x;
  return ((-x2 * r0 + x0 * r2 + t3 - x1 * r2 + x2 * r1 - t6)
      / (-x1 * t8 + t11 + x0 * t8 - t13 + t15 - t16) * t43
      - (t10 * r2 - t10 * r1 - t14 * r2 + r0 * t14 - r0 * t8 + r1 * t8) * t27
          * t33 * x
      + (t11 * r2 - t13 * r1 - t16 * r2 + t6 * t8 + t15 * r0 - t3 * t8) * t27
          * t33);
}

/**
 *  @brief Interpolate a function evaluated at four points using
 *         cubic interpolation
 *
 *  This was generated using Maple's code generator with
 *  optimization.
 */
inline double
interpolate_cubic(double x,
                  double x0, double x1, double x2, double x3,
                  double r0, double r1, double r2, double r3)
{
  double t1;
  double t102;
  double t105;
  double t11;
  double t110;
  double t113;
  double t12;
  double t121;
  double t134;
  double t138;
  double t14;
  double t140;
  double t143;
  double t147;
  double t150;
  double t155;
  double t156;
  double t158;
  double t16;
  double t162;
  double t167;
  double t17;
  double t175;
  double t178;
  double t181;
  double t19;
  double t195;
  double t196;
  double t2;
  double t20;
  double t202;
  double t203;
  double t205;
  double t219;
  double t227;
  double t230;
  double t233;
  double t236;
  double t24;
  double t240;
  double t25;
  double t27;
  double t29;
  double t31;
  double t33;
  double t35;
  double t38;
  double t4;
  double t41;
  double t43;
  double t45;
  double t47;
  double t52;
  double t55;
  double t57;
  double t58;
  double t6;
  double t60;
  double t61;
  double t62;
  double t65;
  double t66;
  double t68;
  double t69;
  double t71;
  double t72;
  double t74;
  double t79;
  double t81;
  double t83;
  double t87;
  double t89;
  double t93;
  double t95;
  double t98;
  t1 = x0 * x0;
  t2 = t1 * x1;
  t4 = t1 * x3;
  t6 = t1 * x2;
  t11 = x0 * r1;
  t12 = x3 * x3;
  t14 = x0 * r2;
  t16 = x1 * x1;
  t17 = x0 * t16;
  t19 = x2 * x2;
  t20 = x0 * t19;
  t24 = -t2 * r3 + t4 * r1 - t6 * r1 - t4 * r2 + t2 * r2 + t6 * r3 - t11 * t12
      + t14 * t12 + t17 * r3 + t20 * r1 - t17 * r2 - t20 * r3;
  t25 = t16 * x3;
  t27 = t16 * x2;
  t29 = r1 * t19;
  t31 = r0 * x2;
  t33 = r0 * t19;
  t35 = r0 * x1;
  t38 = r0 * t16;
  t41 = r1 * x2;
  t43 = r2 * x1;
  t45 = t19 * x1;
  t47 = t25 * r2 - t27 * r3 - t29 * x3 - t31 * t12 + t33 * x3 + t35 * t12
      - t35 * t19 - t38 * x3 + t38 * x2 + t41 * t12 - t43 * t12 + t45 * r3;
  t52 = x0 * x1;
  t55 = x0 * x2;
  t57 = x0 * t12;
  t58 = x2 * x1;
  t60 = x3 * t12;
  t61 = x1 * t12;
  t62 = x2 * t12;
  t65 = x0 * t1;
  t66 = x1 * t65;
  t68 = x1 * t16;
  t69 = x2 * t68;
  t71 = x2 * t19;
  t72 = t71 * x1;
  t74 = t71 * x3;
  t79 = x2 * t65;
  t81 = x0 * t68;
  t83 = t71 * x0;
  t87 = -t66 * r3 + t69 * r0 + t72 * r3 + t74 * r0 - t69 * r3 + t66 * r2
      - t74 * r1 - t79 * r1 + t81 * r3 + t83 * r1 - t83 * r3 + t79 * r3;
  t89 = x3 * t65;
  t93 = t60 * x2;
  t95 = t60 * x0;
  t98 = t60 * x1;
  t102 = x3 * t68;
  t105 = -t72 * r0 - t89 * r2 + t89 * r1 - t81 * r2 - t93 * r0 + t95 * r2
      - t95 * r1 + t98 * r0 - t98 * r2 + t93 * r1 - t102 * r0 + t102 * r2;
  t110 = t68 * t1;
  t113 = t65 * t16;
  t121 = -t79 * t12 + t79 * t16 + t66 * t12 + t110 * x3 - t81 * t12 - t113 * x3
      - t110 * x2 + t83 * t12 - t72 * t12 + t72 * t1 - t74 * t1 + t74 * t16;
  t134 = -t83 * t16 + t69 * t12 + t89 * t19 - t102 * t19 - t66 * t19 + t81 * t19
      - t95 * t19 + t95 * t16 + t98 * t19 - t98 * t1 + t93 * t1 - t93 * t16;
  t138 = t65 * r1;
  t140 = t65 * r2;
  t143 = t65 * t19;
  t147 = t1 * t60;
  t150 = t1 * t71;
  t155 = t138 * t12 - t140 * t12 - t113 * r3 - t143 * r1 + t113 * r2 + t143 * r3
      - t147 * r1 + t147 * r2 + t150 * r1 - t110 * r2 + t110 * r3 - t150 * r3;
  t156 = t16 * t60;
  t158 = t16 * t71;
  t162 = r0 * t68;
  t167 = t19 * t68;
  t175 = -t156 * r2 + t158 * r3 - r1 * t71 * t12 + t162 * t19 - t162 * t12
      + r0 * t71 * t12 - t167 * r3 + t38 * t60 - t33 * t60 - t38 * t71
      + t29 * t60 + r2 * t68 * t12;
  t178 = 0.1e1 / (x0 - x1);
  t181 = t19 * x3;
  t195 = t60 * t19;
  t196 = t71 * t12;
  t202 = t27 * t1 + t181 * t1 - t45 * t1 - t6 * t12 - t25 * t1 + t2 * t12
      + t52 * t71 - t52 * t60 + t55 * t60 - t83 * x3 + t52 * t181 - t17 * t19
      - t58 * t57 + t17 * t12 - t195 + t196 + t58 * t60 - t72 * x3
      + t16 * t19 * x3 - t27 * t12;
  t203 = 0.1e1 / t202;
  t205 = t68 * t12;
  t219 = t181 * t162 - t31 * t205 - t158 * x3 * r0 + t156 * t31 + t35 * t196
      - t35 * t195 - t138 * t181 + t138 * t62 + t150 * x3 * r1 - t147 * t41
      - t11 * t196 + t11 * t195;
  t227 = r3 * t65;
  t230 = r3 * t1;
  t233 = r3 * x0;
  t236 = t140 * t25 - t140 * t61 - t110 * x3 * r2 + t147 * t43 + t14 * t205
      - t156 * t14 - t227 * t27 + t227 * t45 + t230 * t69 - t230 * t72
      - t233 * t167 + t233 * t158;
  t240 = x * x;
  return ((t24 + t47) / (-t45 + t20 + t2 - t17 + t27 - t6)
      / (-t52 * x3 + t52 * x2 - t55 * x3 + t57 - t58 * x3 - t60 + t61 + t62) * x
      * t240 - (t87 + t105) / (t121 + t134) * t240
      + (t155 + t175) * t178 * t203 * x - (t219 + t236) * t178 * t203);
}

/**
 *  @brief Interpolation helper function
 *  @param xi   Value of independent variable at which to evaluate function
 *  @param x    Abscissa
 *  @param r    Function evaluated at abscissa
 */
inline double
interpolate(double xi,
            detran_utilities::vec_dbl &x,
            detran_utilities::vec_dbl &r)
{
  // Preconditions
  Require(x.size() == r.size());
  Require(x.size() > 0);
  Require(x.size() < 5);

  // Number of abscissa for interpolating
  detran_utilities::size_t n = x.size();

  double ri;
  if (n == 1)
    ri = r[0];
  if (n == 2)
    ri = interpolate_linear(xi, x[0], x[1], r[0], r[1]);
  else if (n == 3)
    ri = interpolate_quadratic(xi, x[0], x[1], x[2], r[0], r[1], r[2]);
  else if (n == 4)
    ri = interpolate_cubic(xi, x[0], x[1], x[2], x[3], r[0], r[1], r[2], r[3]);
  else
    THROW("INTERPOLATION WRONG ORDER");
  return ri;
}


} // end namespace erme_response

#endif // erme_response_INTERPOLATION_HH_

//---------------------------------------------------------------------------//
//              end of file Interpolation.hh
//---------------------------------------------------------------------------//
