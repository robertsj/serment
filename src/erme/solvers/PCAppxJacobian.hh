//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PCAppxJacobian.hh
 * \author Jeremy Roberts
 * \date   10/19/2010
 * \brief  An approximate Jacobian matrix for preconditioning.
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 168                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-11-08 14:23:24 -0500 (Tue, 08 Nov 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef PCAPPXJACOBIAN_HH
#define PCAPPXJACOBIANL_HH
#include "LinAlg.hh"
#include "GlobalProblem.hh"

//---------------------------------------------------------------------------//
/*!
 * \class PCAppxJacobian
 * \brief An approximate Jacobian for preconditioning.
 *
 * For effective JFNK, a preconditioner is needed for the linear solves. A
 * straightforward and effective approach is to use an approximate Jacobian.
 * The approach here is to construct a partial Jacobian from the initial
 * guess, and to use it (or its factorization) for the rest of the problem.
 *
 * The approximate Jacobian is defined as
 *  \f[
 *  \mathbf{f'(x)} = \left [\begin{array}{ccc}
 *    (\mathbf{M}\mathbf{R}-\lambda \mathbf{I})                      &
 *          0                                                         &
 *                \mathbf{J_-}                                         \\
 *    (\mathbf{F}-k\mathbf{L})                                       &
 *          -\mathbf{L} \mathbf{J_-}                                  &
 *                0                                                    \\
 *    \mathbf{J^T_-}                                                 &
 *          0                                                         &
 *                0
 *  \end{array}
 *  \right ]  \, ,
 *  \label{eq:jacobian}
 *  \f]
 * where \f$ \mathbf{J_-} \f$, \f$ k \f$, and \f$ \lambda \f$ are the initial
 * guess, likely found from a single crude power iteration.
 *
 * Notice this approximate Jacobian is almost complete, lacking only the
 * finite differences for the derivatives with respect to \f$ k \f$.  No
 * studies have been performed to assess how much dropping these terms affects
 * the convergence of linear solves (or how much time is saved by dropping
 * the terms).  Studies have shown this matrix is an effective preconditioner
 * when coupled with incomplete factorization.
 *
 */
class PCAppxJacobian : public SermentMatrixCRS
{
  // inherits:
  //  mat       M
  //  integer   m, n

public:


  /*!
   * \brief Constructs an approximate Jacobian for preconditioning.
   *
   * \param     a
   * \param     b
   * \param     c
   * \param     nzz
   * \param     unk
   * \param     pr
   *
   */
  PCAppxJacobian(integer a,
                 integer b,
                 integer c,
                 integer const nnz[],
                 SermentVector *unk,
                 GlobalProblem *pr);

  /*!
   * \brief Default destructor.
   */
  ~PCAppxJacobian(){}

  // inherits
  //   void MatVec( SermentVector x, SermentVector y )
  //   void MatVec( SermentVector x )   
  SermentVector *x; // unknown vector

  GlobalProblem *problem; // pointer to global problem (and responses, etc.)

  void myMatVec(Mat &M, Vec &f, Vec &fpf);

private:

};

#endif // PCAPPXJACOBIAN_HH

//---------------------------------------------------------------------------//
//                 end of PCAppxJacobian.hh
//---------------------------------------------------------------------------//

