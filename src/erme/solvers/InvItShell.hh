//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   InvItShell.hh
 * \author Jeremy Roberts
 * \date   10/19/2010
 * \brief  Action of MR-lambdaI
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 140                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-09-14 12:53:40 -0400 (Wed, 14 Sep 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef InvItShell_HH
#define InvItShell_HH
#include "petscsnes.h"
#include "LinAlg.hh"
#include "GlobalProblem.hh"
#include "InvItMatVecWrap.hh"

//===========================================================================//
/*!
 * \class InvItShell
 * \brief 
 *
 */
//===========================================================================//

class InvItShell : public SermentMatrixShell
{
  // inherits:
  //  mat       M
  //  integer   m, n

  public:
    InvItShell( integer a, integer b, void *ctx, GlobalProblem *pr ) 
      : SermentMatrixShell(a,b,this), problem(pr)
    {
        k=0;  
        lambda=0;
        // must set the  matop
        MatShellSetOperation( M, MATOP_MULT, (void(*)())InvItMatVecWrap );
    }; // constructor
   ~InvItShell(){};   // default destructor
    scalar k;
    scalar lambda;
    // inherits
    //   void MatVec( SermentVector x, SermentVector y )
    //   void MatVec( SermentVector x ) 
    GlobalProblem *problem;  // pointer to global problem (and responses, etc.)
    void myMatVec( Mat &M, Vec &x, Vec &y );
    void updateEigs( scalar kval, scalar lval )
    {
        k      = kval;
        lambda = lval;
    }
                            
};

#endif // InvItShell_HH

//---------------------------------------------------------------------------//
//                 end of InvItShell.hh
//---------------------------------------------------------------------------//

