//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SermentMatrixBCRS.hh
 * \author Jeremy Roberts
 * \date   11/23/2010
 * \brief  A SermentMatrix subclass implementing a block CRS matrix
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 75                                            $:Rev of last commit
// $Author:: bert                                       $:Author of last commit
// $Date:: 2010-12-19 09:47:01 -0500 (Sun, 19 Dec 2010) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef SERMENTMATRIXBCRS_HH
#define SERMENTMATRIXBCRS_HH
#include "petscmat.h"
#include "typedefs.hh"
#include "SermentMatrixFull.hh"
#include "SermentVector.hh"

//===========================================================================//
/*!
 * \class SermentMatrixBCRS
 * \brief This concrete matrix class encapsulates a Petsc block CSR matrix.
 *
 *  SermentMatrixCRS builds off SermentMatrixFull.  It
 *  implements methods specific to a complete matrix construction of a 
 *  compressed row storage matrix.
 *
 */
//===========================================================================//

class SermentMatrixBCRS : public SermentMatrixFull
{
  // inherits:
  //  mat       M
  //  integer   m, n, nz


  public:
    SermentMatrixBCRS( integer a, integer b, integer c, 
                       integer d );                       // constructor
   ~SermentMatrixBCRS();                                  // default destructor
      integer bs; // block size
  // inherits
  //   void MatVec( SermentVector x, SermentVector y )
  //   void MatVec( SermentVector x )            


  //  Add values, virtual, implementation depends on type
  //	mat 	    - the matrix
  //	v 	        - a logically two-dimensional array of values
  //	m, idxm 	- the number of block rows and their global block indices
  //	n, idxn 	- the number of block columns and their global block indices
  //	addv 	    - INSERT_VALUES,
  void insertVals(  scalar      values[], 
                    integer     numrow,
                    integer     idxrow[],
                    integer     numcol,
                    integer     idxcol[] );
            
};

#endif // SERMENTMATRIXBCRS_HH

//---------------------------------------------------------------------------//
//                 end of SermentMatrixBCRS.hh
//---------------------------------------------------------------------------//

