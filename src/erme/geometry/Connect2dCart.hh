//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Connect2dCart.hh
 * \author Jeremy Roberts
 * \date   10/19/2010
 * \brief  Connect2dCart class definition.
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 171                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-11-19 20:36:52 -0500 (Sat, 19 Nov 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef CONNECT2DCART_HH
#define CONNECT2DCART_HH
#include <vector>
#include "petscmat.h"
#include "LinAlg.hh"
#include "ConnectMatrix.hh"
using namespace std;

//===========================================================================//
/*!
 * \class Connect2dCart
 * \brief The connectivity matrix for a 2d Cartesian problem.
 *
 *  Connect2dCart builds off SermentMatrixCRS and is the "connectivity" matrix
 *  for 2d Cartesian problems.  The matrix consists at most of one element per 
 *  row and column, and these elements are either 1 or -1.  A nonzero value 
 *  connects like orders of outgoing and incoming currents of adjacent cells.
 *  The -1's arise for odd orders at a reflective boundary.
 *
 */
//===========================================================================//

class Connect2dCart : public ConnectMatrix
{

  /// Typedefs
  //\{
  typedef typename util::SP<Connect2dCart> SP_M;
  //\}

  public:
//  integer numGroups;                    // number of groups
//  integer expOrder;                     // expansion order
//  integer numElements;                  // number of elements
//  integer boundCond[4];                 // boundary conditions
//  scalar  *littleM;                     // zeroth order connectivity
//  integer *littleMidxx;                 // row index for littleM
//  integer *littleMidxy;                 // col index for littleM
//  integer  littleMsize;                 // size of littleM
//  vector<vector<integer> >  elements;

  public:
    Connect2dCart( GlobalInput::SP_globalinput in );                 // constructor
   ~Connect2dCart();                                  // default destructor
  
  // inherits
  //   void MatVec( SermentVector x, SermentVector y )
  //   void MatVec( SermentVector x )            
  //  void insertVals(  scalar      values[], 
  //                    integer     numrow,
  //                    integer     idxrow[],
  //                    integer     numcol,
  //                    integer     idxcol[] );
    void buildMe();
    integer elSum( integer ii, integer II, integer j1, integer j2 );
};

#endif // CONNECT2DCART_HH

//---------------------------------------------------------------------------//
//                 end of Connect2dCart.hh
//---------------------------------------------------------------------------//

