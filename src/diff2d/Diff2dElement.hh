//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Diff2dElement.hh
 * \author Jeremy Roberts
 * \date   10/26/2010
 * \brief  A class for handling input for diff2d problems.
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 177                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-12-09 16:31:49 -0500 (Fri, 09 Dec 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef DIFF2DELEMENT_HH
#define DIFF2DELEMENT_HH

#include <iostream>
#include <fstream>
#include <vector>
#include "petscvec.h"
#include "petscsys.h"
#include "../linalg/typedefs.hh"
#include "../utilities/LegendrePoly.hh"
using namespace std;

//===========================================================================//
/*!
 * \class Diff2dElement
 * \brief A class for representing the domain of an individual diff2d problem.
 *
 *  Diff2dElement is a class that represents geometry and other domain-specific
 *  aspects of a problem.  A Diff2dProblem consists of one or several such
 *  elements.  This arrangement helps streamline response function generation.
 *
 */
//===========================================================================//

class Diff2dElement
{
  public:
        Diff2dElement(){}; 
        ~Diff2dElement(){};
        void destroy()
        {
            VecDestroy( &dx );
            VecDestroy( &dy );
            VecDestroy( &dv );
            P.destroy();
        }
        /*! \brief my identification number (numbering in order of read)
        */
        integer id;
        /*! \brief my description
        */
        string desc;
        /*! \brief volumetric source by coarse mesh and group [][][]
        */
        vector<vector<vector<scalar> > > src;
        /*! \brief material placement by coarse mesh [][]
        */
        vector<vector<integer> >  mt;
        /*! \brief number of x and y fine and coarse meshes *per element*
        */
        integer  nxcm, nxfm, nycm, nyfm;
        /*! \brief x and y coarse mesh (should be *per element*) []
        */
        vector<scalar>  xcm, ycm;
        /*! \brief x and y fine mesh (counts) *per element* []
        */
        vector<integer>  xfm, yfm;
        /*! \brief dx, dy, and dv values
        */
        Vec  dx, dy, dv;
        /*! \brief left boundary   (0=vacuum, 1=reflect, 2=inc cur)
        */
        integer  bcl;
        /*! \brief right boundary  (0=vacuum, 1=reflect, 2=inc cur)
        */
        integer  bcr;
        /*! \brief bottom boundary (0=vacuum, 1=reflect, 2=inc cur)
        */
        integer  bcb;
        /*! \brief top boundary    (0=vacuum, 1=reflect, 2=inc cur)
        */
        integer  bct;
        /*! \brief symmetry    (0=none, 2=two-way, 180 deg, 4=four way)
        */
        integer  sym;
        /*! LEGENDRE POLYNOMIAL CLASS 
        */
        LegendrePoly P;
};

#endif // DIFF2DELEMENT_HH

//---------------------------------------------------------------------------//
//                 end of Diff2dElement.hh
//---------------------------------------------------------------------------//

