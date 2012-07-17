//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SermentMatrixFull.hh
 * \author Jeremy Roberts
 * \date   10/19/2010
 * \brief  A SermentMatrix subclass for full matrix operations
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 54                                            $:Rev of last commit
// $Author:: bert                                       $:Author of last commit
// $Date:: 2010-11-27 17:06:19 -0500 (Sat, 27 Nov 2010) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef SERMENTMATRIXFULL_HH
#define SERMENTMATRIXFULL_HH
#include "petscmat.h"
#include "typedefs.hh"
#include "SermentMatrix.hh"
#include "SermentVector.hh"

//===========================================================================//
/*!
 * \class SermentMatrixFull
 * \brief An abstract matrix class encapsulating a full (constructed) matrix
 *
 *  SermentMatrixFull is an abstract class that builds off SermentMatrix.  It
 *  specifies additional methods specific to a complete matrix construction, 
 *  e.g. the case of a sequential compressed row storage matrix.
 *
 */
//===========================================================================//

class SermentMatrixFull : public SermentMatrix
{
  // inherits:
  //  mat       M
  //  integer   m, n
  public:
  integer nz; // parameter for number of nonzeros
  bool isReady;

  public:
    SermentMatrixFull( integer a, integer b, integer c ) 
      : SermentMatrix(a,b), nz(c) 
    {
    	isReady=false;
    };
   ~SermentMatrixFull(){}; // default destructor
  
  // inherits
  //   void MatVec( SermentVector x, SermentVector y )
  //   void MatVec( SermentVector x )            
  // which can be fully specified using the Petsc MatMult routine.
  void matVec( SermentVector &x, SermentVector &y );
  void matVec( SermentVector x );

  //  Add values, virtual, implementation depends on type
  //	mat 	    - the matrix
  //	v 	        - a logically two-dimensional array of values
  //	m, idxm 	- the number of block rows and their global block indices
  //	n, idxn 	- the number of block columns and their global block indices
  //	addv 	    - INSERT_VALUES,
  virtual void insertVals( scalar      values[], 
                           integer     numrow,
                           integer     idxrow[],
                           integer     numcol,
                           integer     idxcol[] ) = 0;
                            
  void checkReady();
  void viewMe();
};

#endif // SERMENTMATRIXFULL_HH

//---------------------------------------------------------------------------//
//                 end of SermentMatrixFull.hh
//---------------------------------------------------------------------------//

