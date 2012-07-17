//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SermentMatrixCRS.cc
 * \author Jeremy Roberts
 * \date   10/22/2010
 * \brief  Member definitions of concrete class SermentMatrixCRS
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 81                                            $:Rev of last commit
// $Author:: bert                                       $:Author of last commit
// $Date:: 2011-04-03 21:37:16 -0400 (Sun, 03 Apr 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#include "petscvec.h"
#include "petscmat.h"
#include "SermentMatrixCRS.hh"
#include "SermentVector.hh"
#include "typedefs.hh"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
SermentMatrixCRS::SermentMatrixCRS( integer a, integer b, integer c)
  : SermentMatrixFull(a,b,c)
{
	MatCreateSeqAIJ( PETSC_COMM_SELF, m, n, nz, PETSC_NULL, &M );
}

SermentMatrixCRS::SermentMatrixCRS( integer a, integer b, integer c, 
                                    const integer nnz[] )
  : SermentMatrixFull(a,b,c)
{
	MatCreateSeqAIJ( PETSC_COMM_SELF, m, n, nz, nnz, &M );
}

//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
SermentMatrixCRS::~SermentMatrixCRS()
{
    // nothing here right now, as M is destroyed explicitly in Release()
    return; 
}


//---------------------------------------------------------------------------//
/*!
 * \brief Insert values
 *
 * This function adds values to the matrix given a 2-D array of values, the
 * corresponding row and column counts, and indices for placing those values.
 *
 */
void SermentMatrixCRS::insertVals(  scalar     values[], 
                                    integer    numrow,
                                    integer    idxrow[],
                                    integer    numcol,
                                    integer    idxcol[] )
{
    isReady = false; // changing values = no longer ready
    MatSetValues( M, numrow, idxrow, numcol, idxcol, values, 
                  INSERT_VALUES );
    
}

//---------------------------------------------------------------------------//
/*!
 * \brief Insert value
 *
 * This function adds a single values to the matrix given an index.
 *
 */
void SermentMatrixCRS::insertVal(  scalar      value, 
                                   integer     row,
                                   integer     col )
{
    isReady = false; // changing values = no longer ready
    MatSetValue( M, row, col, value, INSERT_VALUES );
    
}
//---------------------------------------------------------------------------//
//                 end of SermentMatrixCRS.cc
//---------------------------------------------------------------------------//

