//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SermentMatrixCRS.hh
 * \author Jeremy Roberts
 * \date   10/19/2010
 * \brief  A SermentMatrix subclass implementing a CRS matrix
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 68                                            $:Rev of last commit
// $Author:: bert                                       $:Author of last commit
// $Date:: 2010-12-09 17:17:34 -0500 (Thu, 09 Dec 2010) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef SERMENTMATRIXCRS_HH
#define SERMENTMATRIXCRS_HH
#include "petscmat.h"
#include "typedefs.hh"
#include "SermentMatrixFull.hh"
#include "SermentVector.hh"

//===========================================================================//
/*!
 * \class SermentMatrixCRS
 * \brief This concrete matrix class encapsulates a Petsc CSR matrix.
 *
 *  SermentMatrixCRS builds off SermentMatrixFull.  It
 *  implements methods specific to a complete matrix construction of a 
 *  compressed row storage matrix.
 *
 */
//===========================================================================//

class SermentMatrixCRS : public SermentMatrixFull
{
  // inherits:
  //  mat       M
  //  integer   m, n, nz

  public:
    SermentMatrixCRS( integer a, integer b, integer c ); // constructor
    SermentMatrixCRS( integer a, integer b, integer c, 
                      const integer nnz[] );             // constructor2
   ~SermentMatrixCRS();                                  // default destructor
  
  // inherits
  //   void MatVec( SermentVector x, SermentVector y )
  //   void MatVec( SermentVector x )            


  //  Add values
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

  //  Add value
  //	mat 	    - the matrix
  //	v 	        - a logically two-dimensional array of values
  //	m           - the row
  //	n         	- the column
  //	addv 	    - INSERT_VALUES,
  void insertVal(  scalar      value, 
                   integer     row,
                   integer     col );
            
};

#endif // SERMENTMATRIXCRS_HH

//---------------------------------------------------------------------------//
//                 end of SermentMatrixCRS.hh
//---------------------------------------------------------------------------//

