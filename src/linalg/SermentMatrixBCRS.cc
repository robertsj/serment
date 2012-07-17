//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SermentMatrixBCRS.cc
 * \author Jeremy Roberts
 * \date   11/23/2010
 * \brief  Member definitions of concrete class SermentMatrixBCRS
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 75                                            $:Rev of last commit
// $Author:: bert                                       $:Author of last commit
// $Date:: 2010-12-19 09:47:01 -0500 (Sun, 19 Dec 2010) $:Date of last commit
//---------------------------------------------------------------------------//

#include "petscvec.h"
#include "petscmat.h"
#include "SermentMatrixBCRS.hh"
#include "SermentVector.hh"
#include "typedefs.hh"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
SermentMatrixBCRS::SermentMatrixBCRS( integer a, integer b, integer c, 
                                      integer d )
  : SermentMatrixFull(a,b,c), bs(d)
{
	MatCreateSeqBAIJ( PETSC_COMM_SELF, bs, m, n, nz, PETSC_NULL, &M );
    // bs = block size (i.e. length of one row)
    // m  = number of (block) rows
    // n  = number of (block) cols
    // nz = number of nonzero blocks per row (same for all rows)
    // Using column-oriented construction of the blocks (as apparently, that's
    // how I think best...
    MatSetOption( M, MAT_ROW_ORIENTED, PETSC_FALSE );
}

//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
SermentMatrixBCRS::~SermentMatrixBCRS()
{
    // nothing here right now, as M is destroyed explicitly in Release()
    return; 
}

//---------------------------------------------------------------------------//
/*!
 * \brief This function insert values into a block CSR matrix.
 *
 * This function adds values to the matrix given a 2-D array of values, the
 * corresponding row and column counts, and indices for placing those values.
 * Note, the values are placed in block format, but they themselves are still
 * in a 1d array.  For example, if m=n=2, and we want the following blocks
 *   1  2  | 3  4
 *   5  6  | 7  8
 *   - - - | - - -
 *   9  10 | 11 12
 *   13 14 | 15 16
 * the vector v would have {1, 5, ..., 16} sequentially, i.e. in column-major 
 * form.  Likely, the BCSR matrices to be computed in Serment, namely the
 * response matrix, will be constructed one block at a time, since the current
 * response blocks are stored in separate 1d arrays for each element.
 *
 */
void SermentMatrixBCRS::insertVals(  scalar     values[], 
                                     integer    numrow,
                                     integer    idxrow[],
                                     integer    numcol,
                                     integer    idxcol[] )
{
    isReady = false; // changing values = no longer ready
    MatSetValuesBlocked( M, numrow, idxrow, numcol, idxcol, values, 
                         INSERT_VALUES );
    
}

//---------------------------------------------------------------------------//
//                 end of SermentMatrixBCRS.cc
//---------------------------------------------------------------------------//

