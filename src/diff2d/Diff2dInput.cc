//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Diff2dInput.cc
 * \author Jeremy Roberts
 * \date   10/26/2010
 * \brief  Member definitions of class Diff2dInput
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 87                                            $:Rev of last commit
// $Author:: bert                                       $:Author of last commit
// $Date:: 2011-04-11 13:51:25 -0400 (Mon, 11 Apr 2011) $:Date of last commit
//---------------------------------------------------------------------------//
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <math.h>
#include "petscvec.h"
#include "petscmat.h"
#include "Diff2dInput.hh"
#include "../linalg/typedefs.hh"

using namespace std;

//---------------------------------------------------------------------------//
/*!
 * \brief readInput
 *
 * This function reads the local problem input and initializes all data
 * subsequently used to make the operater matrices and solve the problem.
 * 
 * The input file structure, and consequently, the reader, should be made
 * more general in future refinements.  However, the basic data is all that
 * is needed for now.
 *
 */
void Diff2dInput::readInput( char *file )
{
    debug = false;

    // announce yourself
    cout << "|-----------------------------------------------|" << endl;
    cout << "|                    diff2d                     |" << endl;
    cout << "|           a simple 2-d diffusion code         |" << endl;
    cout << "|                                               |" << endl;
    cout << "|  Version 0.1, 04/09/2011                      |" << endl;
    cout << "|-----------------------------------------------|" << endl;
    cout << " opening user file: " << file << endl << endl;

	//-----------------------------------------------------------------------
    // open file and exit on error
    ifstream inFile( file );
    if (! inFile )
    {
        cout << "ERROR: The file cannot be opened" << endl;
        exit(EXIT_FAILURE);
    }

    string tmp;              // temporary string for reading lines
    istringstream tmpcin;    // stream for outputing tmp into variables
    string tmp2;             // string for storing extras, e.g. "numg="
    skipStuff( inFile );     // skip header stuff


	//-----------------------------------------------------------------------
    // reading GLOBAL CONTROL
    if (debug) cout << "...reading control block..." << endl;
    getline(inFile,tmp); tmpcin.str(tmp); tmpcin >> tmp2; tmpcin >> ptype;
    getline(inFile,tmp); tmpcin.str(tmp); tmpcin >> tmp2; tmpcin >> epsk;
    getline(inFile,tmp); tmpcin.str(tmp); tmpcin >> tmp2; tmpcin >> epss;
    getline(inFile,tmp); tmpcin.str(tmp); tmpcin >> tmp2; tmpcin >> maxit;
    getline(inFile,tmp); tmpcin.str(tmp); tmpcin >> tmp2; tmpcin >> numel;
    if ( ptype == 2 )
        getline(inFile,tmp); tmpcin.str(tmp); tmpcin >> tmp2; tmpcin >> maxOrder;

    if (debug)
    {
    cout << " control parameters: " << endl;
    cout << " ptype=" << ptype << " epsk=" << epsk << " epss=" << epss 
         << " maxit=" << maxit << " numel=" << numel << endl << endl;
    }

    if ( ptype == 2) cout << " maxOrder= " << maxOrder << endl;
    skipStuff( inFile ); 

	//-----------------------------------------------------------------------
    // reading OUTPUT CONTROL
    if (debug) cout << "...reading control block..." << endl;
    getline(inFile,tmp); tmpcin.str(tmp); tmpcin >> tmp2; tmpcin >> printout;
    getline(inFile,tmp); tmpcin.str(tmp); tmpcin >> tmp2; tmpcin >> plotout;
    getline(inFile,tmp); tmpcin.str(tmp); tmpcin >> tmp2; tmpcin >> outfile;
    getline(inFile,tmp); tmpcin.str(tmp); tmpcin >> tmp2; tmpcin >> fluxfile;
    if (debug)
    {
        cout << " control parameters: " << endl;
        cout << " printout=" << printout << " plotout=" << plotout 
             << " outfile=" << outfile << " fluxfile=" << fluxfile << endl << endl;
    }
    skipStuff( inFile );


	//-----------------------------------------------------------------------
    // reading MATERIAL DATA
    //cout << "...reading material block..." << endl;
    getline(inFile,tmp); tmpcin.str(tmp); tmpcin >> tmp2; tmpcin >> numg;
    getline(inFile,tmp); tmpcin.str(tmp); tmpcin >> tmp2; tmpcin >> numm;
    if (debug)
    {
    cout << " numg=" << numg << " numm=" << numm << endl;
    }
    skipStuff( inFile ); 
    vector<vector<scalar> > data; // data [numm*numg][4+numg] 
    data.resize(numm*numg);       // temp array for reading cross-section data
    for(int i = 0; i < numm*numg; ++i)
        data[i].resize(4+numg);

    dc.resize(numm);
    sr.resize(numm);
    ab.resize(numm);
    ns.resize(numm);
    xi.resize(numm);
    sc.resize(numm);
    for(int i = 0; i < numm; ++i)
    {
        dc[i].resize(numg);
        sr[i].resize(numg);
        ab[i].resize(numg);
        ns[i].resize(numg);
        xi[i].resize(numg);
        sc[i].resize(numg);
        for(int j = 0; j < numg; ++j)
            sc[i][j].resize(numg);
    }
	// read data
    for(int i = 0; i < numm*numg; ++i)
    {
        for(int j = 0; j < numg+4; ++j)
			inFile >> data[i][j];
    }
	// verify
    if (debug)
        {
        cout << " ... data = ";
        for(int i = 0; i < numm*numg; ++i)
        {
            cout << endl;
            for(int j = 0; j < numg+4; ++j)
            {
			    cout.width(6); cout << data[i][j] << " ";
            }
        }
	    cout << endl << "-----" << endl;
    }

    skipStuff( inFile );

	//-----------------------------------------------------------------------
    // reading READ THE ELEMENTS
    elements.resize(numel);
    for ( int e = 0; e < numel; ++e )
    {
        if (debug) cout << " element: " << e << endl;
        elements[e].id = e; // each element gets a number assigned to it
        skipStuff( inFile );
        getline(inFile,tmp); tmpcin.str(tmp); tmpcin >> tmp2; tmpcin >> elements[e].desc;
        getline(inFile,tmp); tmpcin.str(tmp); tmpcin >> tmp2; tmpcin >> elements[e].nxcm;
        getline(inFile,tmp); tmpcin.str(tmp); tmpcin >> tmp2; tmpcin >> elements[e].nxfm;
        getline(inFile,tmp); tmpcin.str(tmp); tmpcin >> tmp2; tmpcin >> elements[e].nycm;
        getline(inFile,tmp); tmpcin.str(tmp); tmpcin >> tmp2; tmpcin >> elements[e].nyfm;
        getline(inFile,tmp); tmpcin.str(tmp); tmpcin >> tmp2; tmpcin >> elements[e].bcl;
        getline(inFile,tmp); tmpcin.str(tmp); tmpcin >> tmp2; tmpcin >> elements[e].bcr;
        getline(inFile,tmp); tmpcin.str(tmp); tmpcin >> tmp2; tmpcin >> elements[e].bcb;
        getline(inFile,tmp); tmpcin.str(tmp); tmpcin >> tmp2; tmpcin >> elements[e].bct;
        getline(inFile,tmp); tmpcin.str(tmp); tmpcin >> tmp2; tmpcin >> elements[e].sym;
        if (debug)
        {
            cout << " basic parameters: " << endl;
            cout << " desc= " << elements[e].desc << " id= " << elements[e].id << endl;
            cout << " numg=" << numg << " numm=" << numm << " nxcm=" << elements[e].nxcm 
                 << " nxfm=" << elements[e].nxfm << " nycm=" << elements[e].nycm << " nyfm=" << elements[e].nyfm << endl;
            cout << " bcl=" << elements[e].bcl << " bcr=" << elements[e].bcr 
                 << " bcb=" << elements[e].bcb << " bct=" << elements[e].bct << " sym=" << elements[e].sym << endl << endl;
        }

        skipStuff( inFile );

        // allocate class member in order to read array values
        elements[e].xcm.resize(elements[e].nxcm+1);
        elements[e].ycm.resize(elements[e].nycm+1);
        elements[e].xfm.resize(elements[e].nxcm);
        elements[e].yfm.resize(elements[e].nycm);
	    // src [nxcm][nycm]
        elements[e].src.resize(elements[e].nxcm);
        for(int i = 0; i < elements[e].nxcm; ++i)
        {
            elements[e].src[i].resize(elements[e].nycm);
            for(int j = 0; j < elements[e].nycm; ++j)
                elements[e].src[i][j].resize(numg);
        }
        // mt [nxcm][nycm]
        elements[e].mt.resize(elements[e].nxcm);
        for(int i = 0; i < elements[e].nxcm; ++i)
            elements[e].mt[i].resize(elements[e].nycm);

	    // read geometry and material placement

	    // x coarse mesh
        for(int i = 0; i < elements[e].nxcm+1; ++i)
		    inFile >> elements[e].xcm[i];
        if (debug)
        {
            cout << "...xcm = " << endl;
            for(int i = 0; i < elements[e].nxcm+1; ++i)
		        cout << elements[e].xcm[i] << " ";
            cout << endl << endl;
        }
	    skipStuff( inFile );

	    // x fine mesh
        for(int i = 0; i < elements[e].nxcm; ++i)
		    inFile >> elements[e].xfm[i];
        if (debug)
        {
            cout << "...xfm = " << endl;
            for(int i = 0; i < elements[e].nxcm; ++i)
		        cout << elements[e].xfm[i] << " ";
            cout << endl << endl;
        }
	    skipStuff( inFile );

	    // y coarse mesh
        for(int i = 0; i < elements[e].nycm+1; ++i)
		    inFile >> elements[e].ycm[i];
        if (debug)
        {  
            cout << "...ycm = " << endl;
            for(int i = 0; i < elements[e].nycm+1; ++i)
		        cout << elements[e].ycm[i] << " ";
            cout << endl << endl;
        }
	    skipStuff( inFile );

	    // y fine mesh
        for(int i = 0; i < elements[e].nycm; ++i)
		    inFile >> elements[e].yfm[i];

        if (debug)
        {
            cout << "...yfm = " << endl;
            for(int i = 0; i < elements[e].nycm; ++i)
		        cout << elements[e].yfm[i] << " ";
            cout << endl << endl;
        }

	    skipStuff( inFile );

	    // read material block
        for(int i = 0; i < elements[e].nxcm; ++i)
        {
            for(int j = 0; j < elements[e].nycm; ++j)
		    {
			    inFile >> elements[e].mt[i][j];
                --elements[e].mt[i][j];
		    }
        }

        // read volume source *if* needed
        if ( ptype < 1 )
        {
            for(int g = 0; g < numg; ++g)
            {
	            skipStuff( inFile );
                for(int i = 0; i < elements[e].nxcm; ++i)
                {
                    for(int j = 0; j < elements[e].nycm; ++j)
			            inFile >> elements[e].src[i][j][g];
                }
            }
        }
    } // end loop over elements

	//-----------------------------------------------------------------------
    //  produce the data vectors for later matrix construction
    for(int m = 0; m < numm; ++m)
    {
        for(int g = 0; g < numg; ++g)
		{
            cout << " m = " << m << " g = " << g << endl;
            dc[m][g] = data[numg*m+g][0];  // diffusion coefficient
            sr[m][g] = data[numg*m+g][1];  // removal cross-section
            ns[m][g] = data[numg*m+g][2];  // fission cross-section
            xi[m][g] = data[numg*m+g][3];  // fission spectrum

        	for(int gg = 0; gg < numg; ++gg)
				sc[m][g][gg] = data[numg*m+g][4+gg];
            if ( numg >= 1 )
			{
				scalar sum = 0;
				for(int i = numg*m; i < numg*(m+1); ++i)
					sum = sum + data[i][4+g];
            	ab[m][g] = sr[m][g] - sum;
			}
            else 
			{
            	ab[m][g] = sr[m][g];
			}
		}
    }

    if (debug)
    {
        // verify mt and src
        for(int e = 0; e < numel; ++e)
        {    cout << "...reading material block..." << endl;
            for(int m = 0; m < elements[e].nxcm; ++m)
            {
                for(int g = 0; g < elements[e].nycm; ++g)
                    cout << elements[e].mt[m][g] << " ";
                cout << " --- " << endl;
            }
            if ( ptype < 1 )
            {
                for(int g = 0; g < numg; ++g)
                {
	                cout << " source group " << g+1 << endl;
                    for(int i = 0; i < elements[e].nxcm; ++i)
                    {
                        for(int j = 0; j < elements[e].nycm; ++j)
			                cout <<  elements[e].src[i][j][g] << " ";
                        cout << endl;
                    }
                }
            }
        }

        // verify materials
        for ( int m = 0; m < numm; m++ )
        {
            cout << "       material: " << m << endl;
            for ( int g = 0; g < numg; g++ )
            {
                cout << "          group: " << g << endl;
                cout << "                   ";
                printf (" %10.6f ", dc[m][g] );
                printf (" %10.6f ", sr[m][g] );
                printf (" %10.6f ", ab[m][g] );
                printf (" %10.6f ", ns[m][g] );
                printf (" %10.6f ", xi[m][g] );
                for ( int gg = 0; gg < numg; gg++ )
                {
                    printf (" %10.6f ", sc[m][g][gg] );
                }
                cout << endl;

            }
            cout << endl;
        }
    }

    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief This function skips over blocks of comments and blank lines.
 *
 */
void Diff2dInput::skipStuff( ifstream &in )
{
    // get past comment blocks and spaces
    string tmp;
    int n = 0;
    while( in.peek() == '%' or in.peek() == '\n' )
    {        
        getline(in,tmp);
        n++;
    }
    if (debug) cout << n << " lines skipped " << endl << endl;
    return;
}

//---------------------------------------------------------------------------//
//                 end of Diff2dInput.cc
//---------------------------------------------------------------------------//

