//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SermentMatrixShell.hh
 * \author Jeremy Roberts
 * \date   10/19/2010
 * \brief  A SermentMatrix subclass for shell matrix operations
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 63                                            $:Rev of last commit
// $Author:: bert                                       $:Author of last commit
// $Date:: 2010-12-06 19:34:00 -0500 (Mon, 06 Dec 2010) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef SERMENTMATRIXSHELL_HH
#define SERMENTMATRIXSHELL_HH
#include "petscmat.h"
#include "typedefs.hh"
#include "SermentMatrix.hh"
#include "SermentVector.hh"

//===========================================================================//
/*!
 * \class SermentMatrixShell
 * \brief This is an abstract class encapsulating a shell matrix (action only).
 *
 *  SermentMatrixShell is an abstract class that builds off SermentMatrix.  It
 *  specifies additional methods specific for the internal mechanism of
 *  matrix-vector operations
 *
 */
//===========================================================================//

class SermentMatrixShell : public SermentMatrix
{
  // inherits:
  //  mat       M
  //  integer   m, n

  public:
    SermentMatrixShell( integer a, integer b, void *ctx );
   ~SermentMatrixShell(){};   // default destructor
  
  // inherits
  //   void MatVec( SermentVector x, SermentVector y )
  //   void MatVec( SermentVector x )            
  // which can be fully specified using the internal MyMatVec routine
  void matVec( SermentVector &x, SermentVector &y );
  void matVec( SermentVector x ){return;}; // not yet implemented
  virtual void myMatVec( Mat &A, Vec &x, Vec &y ) = 0;
                            
};

#endif // SERMENTMATRIXSHELL_HH

//---------------------------------------------------------------------------//
//                 end of SermentMatrixShell.hh
//---------------------------------------------------------------------------//

