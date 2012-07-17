//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Diff2dOutput.cc
 * \author Jeremy Roberts
 * \date   10/26/2010
 * \brief  Member definitions of class Diff2dOutput
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 177                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-12-09 16:31:49 -0500 (Fri, 09 Dec 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"

#include "Diff2dOutput.hh"


using namespace std;

#ifdef SERMENT_ENABLE_SILO
//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
Diff2dOutput::Diff2dOutput( Diff2dInput &in )
{
    inp = &in;
    silocreated = false;
    file = NULL;
    return;
}

//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
Diff2dOutput::~Diff2dOutput()
{

    if (silocreated)
    {
        cout << " closing silo file " << endl;
        DBClose(file);
    }
    cout << " --- killing the Output " << endl;
    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Close the silo filo---used for verification
 *
 */
void Diff2dOutput::closeSilo()
{
    if (silocreated)
    {
        cout << " closing silo file " << endl;
        DBClose(file);
        silocreated = false;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Control output
 *
 */
void Diff2dOutput::doOutput( Diff2dProblem &prob, Diff2dSolver &sol )
{
    
    if ( inp->printout == 1 )
    {
        if ( inp->ptype == 2 ) cout << " warning: ptype=2 w/ standard output " << endl;
        printFlux();
    }
    if ( inp->plotout == 1 )
    {
        if ( inp->ptype == 2 ) cout << " warning: ptype=2 w/ silo! " << endl;
        plotFlux( prob, sol );
    }
    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Print the flux and other things to ASCII
 *
 */
void Diff2dOutput::printFlux()
{
    cout << " printing output ... not really. " << endl;
    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Create a silo file for plotting in Visit
 *
 */
void Diff2dOutput::plotFlux( Diff2dProblem &prob, Diff2dSolver &sol )
{
    cout << " PLOTTING FLUX " << endl;
    

    char *coordnames[2];        /* Names of the coordinates */        
    float *coords[2];           /* The array of coordinate arrays */
    int dims[2];                /* The number of nodes in each dimension */

    // create the silo file; *one* file for all elements, but different
    // elements will have different meshes.
    if ( !silocreated )
    {
        cout << " creating silo " << endl;
        file = DBCreate( (char *)inp->fluxfile.c_str(), DB_CLOBBER, 
                         DB_LOCAL, "group fluxes", DB_PDB );
        silocreated = true;
        if ( !file )
        {
            cout << " ERROR: can't create silo file! " << endl;
            return;
        }
    }
    // name the coordinates
    char c1[5] = "x"; coordnames[0] = c1;
    char c2[5] = "y"; coordnames[1] = c2;
    // set coordinate dimensions
    dims[0] = prob.el->nxfm+1;
    dims[1] = prob.el->nyfm+1;
    // create the coordinate vectors (mesh edge!)
    float x[dims[0]], y[dims[1]];
    x[0]=prob.el->xcm[0];
    y[0]=prob.el->ycm[0];
    for ( int i = 1; i < dims[0]; i++ )
        x[i]=x[i-1]+prob.dx[i-1];
    for ( int i = 1; i < dims[1]; i++ )
        y[i]=y[i-1]+prob.dy[i-1]; 
    coords[0] = x;
    coords[1] = y;

    // make the mesh itself for the specific element
    char buffer1 [12];
    sprintf(buffer1,"mesh_elem%i",prob.el->id);
    cout << " BUFFER1 = " << buffer1 << endl;
    DBPutQuadmesh( file, buffer1, NULL, coords, dims, 2, 
                   DB_FLOAT, DB_COLLINEAR, NULL);

    scalar *tmp;
    char buffer2 [14];
    dims[0]--;
    dims[1]--;
    // get the scalar values from the fluxes
    for (integer g = 0; g < inp->numg; g++)
    {
        sprintf(buffer2,"elem%igroup%i",prob.el->id,g+1);
        // get the flux
        VecGetArray( sol.phi[g], &tmp );
        DBPutQuadvar1( file, buffer2, buffer1, tmp, dims, 2, NULL, 0, 
                       DB_DOUBLE, DB_ZONECENT, NULL);
        // restore the flux
        VecRestoreArray( sol.phi[g], &tmp );
    }

    return;
}
#else
// temporary output stubs for use without silo... \todo create a better output.
Diff2dOutput::Diff2dOutput( Diff2dInput &in )
{
    return;
}
Diff2dOutput::~Diff2dOutput()
{
  return;
}
void Diff2dOutput::closeSilo()
{
  return;
}
void Diff2dOutput::doOutput( Diff2dProblem &prob, Diff2dSolver &sol )
{
  return;
}
void Diff2dOutput::printFlux()
{
  return;
}
void Diff2dOutput::plotFlux( Diff2dProblem &prob, Diff2dSolver &sol )
{
  return;
}
#endif
//---------------------------------------------------------------------------//
//                 end of Diff2dOutput.cc
//---------------------------------------------------------------------------//

