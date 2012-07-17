//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   JacobianShell.hh
 * \author Jeremy Roberts
 * \date   10/19/2010
 * \brief  Action of the Jacobian.
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 140                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-09-14 12:53:40 -0400 (Wed, 14 Sep 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef JACOBIANSHELL_HH
#define JACOBIANSHELL_HH
#include "petscsnes.h"
#include "LinAlg.hh"
#include "GlobalProblem.hh"
#include "MyMatVecWrap.hh"

//===========================================================================//
/*!
 * \class JacobianShell
 * \brief 
 *
 */
//===========================================================================//

class JacobianShell : public SermentMatrixShell
{
  // inherits:
  //  mat       M
  //  integer   m, n

  public:
    JacobianShell( integer a, integer b, void *ctx, SermentVector *unk, 
                    SermentVector *res, GlobalProblem *pr, SNES &sn ) 
      : SermentMatrixShell(a,b,ctx), x(unk), ff(res), problem(pr), snes(sn) 
    {
        // must set the  matop
        MatShellSetOperation( M, MATOP_MULT, (void(*)())myMatVecWrap );
    }; // constructor
   ~JacobianShell(){};   // default destructor
  
  // inherits
  //   void MatVec( SermentVector x, SermentVector y )
  //   void MatVec( SermentVector x ) 
  SermentVector *x;        // unknown vector          
  SermentVector *ff;        // pointer to residual      
  GlobalProblem *problem;  // pointer to global problem (and responses, etc.)
  SNES snes;
  void myMatVec( Mat &M, Vec &f, Vec &fpf );
  
  // *temporary* fix for jacobian rf-updates
  // for each new k, update the k1 and k2 R, L, F, and A once and store
  // a personal copy

                          
};

#endif // JACOBIANSHELL_HH

//---------------------------------------------------------------------------//
//                 end of JacobianShell.hh
//---------------------------------------------------------------------------//

