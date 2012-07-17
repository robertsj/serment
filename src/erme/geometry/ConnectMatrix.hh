//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ConnectMatrix.hh
 * \author Jeremy Roberts
 * \date   10/19/2010
 * \brief  ConnectMatrix class definition.
 * \note   Copyright (C) 2011 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 171                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-11-19 20:36:52 -0500 (Sat, 19 Nov 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef ConnectMatrix_HH
#define ConnectMatrix_HH
#include <vector>
#include "petscmat.h"
#include "LinAlg.hh"
#include "GlobalInput.hh"
using namespace std;

//===========================================================================//
/*!
 * \class ConnectMatrix
 * \brief The base connectivity matrix.
 *
 *  ConnectMatrix builds off SermentMatrixCRS.
 *
 */
//===========================================================================//

class ConnectMatrix : public SermentMatrixCRS
{
    // inherits:
    //   mat       M
    //   integer   m, n, nz
    //   void MatVec( SermentVector x, SermentVector y )
    //   void MatVec( SermentVector x )            
    //   void insertVals(  scalar      values[], 
    //                     integer     numrow,
    //                     integer     idxrow[],
    //                     integer     numcol,
    //                     integer     idxcol[] );
    // note: SermentMatrixCRS( integer a, integer b, integer c )

  public:
    virtual void buildMe() = 0;
    integer *getMindex(){ return littleMidx; };
    integer getMsize(){ return littleMsize; };
    ConnectMatrix( GlobalInput *in ) 
        : SermentMatrixCRS( in->degfree, in->degfree, 1 ),
          input(in){}                                      // constructor
   ~ConnectMatrix(){};                                     // default destructor

  protected:
    GlobalInput *input;                   // pointer to input;
    integer *littleMidx;                  // diagonal index for littleM
    integer  littleMsize;                 // size of littleM

//    integer elSum( integer ii, integer II, integer j1, integer j2 );
};

#endif // ConnectMatrix_HH

//---------------------------------------------------------------------------//
//                 end of ConnectMatrix.hh
//---------------------------------------------------------------------------//

