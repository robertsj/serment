//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LegendrePoly.hh
 * \author Jeremy Roberts
 * \date   10/26/2010
 * \brief  A class for handling input for diff2d problems.
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 106                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-06-15 20:35:53 -0400 (Wed, 15 Jun 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef LEGENDREPOLY_HH
#define LEGENDREPOLY_HH

#include <iostream>
#include <fstream>
#include <vector>
#include "petscvec.h"
#include "petscsys.h"
#include "../linalg/typedefs.hh"
using namespace std;

//===========================================================================//
/*!
 * \class LegendrePoly
 * \brief This class is a complete Legendre polynomial basis set for a 
 *        given element.
 *
 *  LegendrePoly computes and stores the (discrete) Legendre polynomials
 *  needed to produce incident boundary conditions and to expand outgoing
 *  partial currents as response functions.  The polynomial vectors are 
 *  produced once at initialization and stored for subsequent application.
 *
 */
//===========================================================================//

class LegendrePoly
{

  public: // probably want these private...public for testing

        /// polynomials for expansions along horizontal edge
        vector<Vec> Px;
        /// polynomials for expansions along vertical edge
        vector<Vec> Py;     
        /// size of Px for a given order
        integer numx;
        /// size of Py for a given order
        integer numy;
        /// maximum order of the expansion (same for X and Y!)
        integer maxOrder;  
        /// weights for horizontal expansion
        vector<scalar> wx;  
        /// weights for vertical expansion
        vector<scalar> wy;
        
  public:
        LegendrePoly();
        ~LegendrePoly();
        void buildMe( integer nx, integer ny, integer maxOrder );
        void makeWeights();
        void expandCur( Vec &partialCurrent, scalar legCoefs[] );
        scalar factorial( integer x );
        void destroy();
//        Vec H( integer order );
//        Vec V( integer order );

};

#endif // LEGENDREPOLY_HH

//---------------------------------------------------------------------------//
//                 end of LegendrePoly.hh
//---------------------------------------------------------------------------//

