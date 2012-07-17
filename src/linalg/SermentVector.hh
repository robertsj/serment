//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SermentVector.hh
 * \author Jeremy Roberts
 * \date   10/15/2010
 * \brief  A Vector class that encapsulates an external (PETSc) vector class.
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 167                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-11-07 23:01:04 -0500 (Mon, 07 Nov 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef SERMENTVECTOR_HH
#define SERMENTVECTOR_HH
#include "petscvec.h"
#include "typedefs.hh"

#include "utilities/SP.hh"

//===========================================================================//
/*!
 * \class SermentVector
 * \brief A Vector class that encapsulates an external (PETSc) vector class.
 *
 *  SermentVector is a class used for all vector quantities that are operated
 *  on by matrices (i.e. SermentMatrices).  SermentVector encapsulates an 
 *  external package's vector to hide external-specific code within the main
 *  body of Serment.  Currently, PETSc is used for all linear algebra.  Also,
 *  only a sequential code is being implemented currently.  However, because
 *  all vector (and matrix) construction is hidded in these encapsulations, a
 *  move to parallel will not change the main code.
 *
 */
//===========================================================================//

class SermentVector
{
public:

  /// Typedefs
  //\{
  typedef typename util::SP<SermentVector> SP_vector;
  //\}

  Vec V; // PETSc vector
  integer m; // number of rows

private:
  bool isReady; // is it ready for operation?

public:
  SermentVector(integer m);
  ~SermentVector();

  void insertVal(integer row, scalar value);
  void insertVals(integer ni, const integer ix[], const scalar y[]);
  void vecSet(scalar a);
  void vecAVPY(scalar a, SermentVector &Y);
  void vecAYPV(scalar a, SermentVector &Y);
  scalar vecDot(SermentVector &Y);
  void vecPointMult(SermentVector &Y, SermentVector &Z);
  void vecScale(scalar a);
  void vecCopy(SermentVector &Y);



  integer Length();
  void releaseMe();
  void checkReady();
  void viewMe();


};

#endif // SERMENTVECTOR_HH
//---------------------------------------------------------------------------//
//                 end of SermentVector.hh
//---------------------------------------------------------------------------//

