//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SermentMatrix.hh
 * \author Jeremy Roberts
 * \date   10/19/2010
 * \brief  A base Matrix class that encapsulates an external matrix
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 177                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-12-09 16:31:49 -0500 (Fri, 09 Dec 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef SERMENTMATRIX_HH
#define SERMENTMATRIX_HH
#include "petscmat.h"
#include "typedefs.hh"
#include "SermentVector.hh"
//===========================================================================//
/*!
 * \class SermentMatrix
 * \brief An abstract matrix class that encapsulates an external matrix class.
 *
 *  SermentMatrix is the base class used for all matrix quantities used in 
 *  SERMENT's linear operations.  The base (i.e. abstract) class and its 
 *  subclasses support one common contract, namely to provide action of a
 *  vector.  Since the ultimate implementation depends on the underlying matrix
 *  format, this is defined as a virtual method with appropriate signatures.
 *
 */
//===========================================================================//

class SermentMatrix
{
  public:
    Mat         M;  // Petsc matrix
    integer     m;  // number of rows
    integer     n;  // number of columns
  public:
    SermentMatrix( integer a, integer b ) : m(a), n(b) {};
   ~SermentMatrix(){};
	virtual void matVec( SermentVector &x, SermentVector &y ) = 0; // Mx --> y
    virtual void matVec( SermentVector x ) = 0;                    // Mx --> x
//---------------------------------------------------------------------------//
/*!
 * \brief Deallocate the memory for the Petsc matrix explicitly.
 */
    void releaseMe(){MatDestroy( &M );};
};

#endif // SERMENTMATRIX_HH

//---------------------------------------------------------------------------//
//                 end of SermentMatrix.hh
//---------------------------------------------------------------------------//

