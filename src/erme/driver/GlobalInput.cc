//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   GlobalInput.cc
 * \author Jeremy Roberts
 * \date   10/26/2010
 * \brief  Member definitions of class GlobalInput
 * \note   Copyright (C) 2011 Jeremy Roberts.  Several functions are based on
           code by K. Huff from the cyclus project.
 */
//---------------------------------------------------------------------------//
// $Rev:: 140                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-09-14 12:53:40 -0400 (Wed, 14 Sep 2011) $:Date of last commit
//---------------------------------------------------------------------------//
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <math.h>
#include <ctime>
#include "petscvec.h"
#include "petscmat.h"
#include "LinAlg.hh"
#include "GlobalInput.hh"

using namespace std;

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
GlobalInput::GlobalInput(string S)
 : InputXML(S)
{
	// defaults for optional variables
    rfsource = 0;
    rfsink   = 0;
    dbfile   = "n/a";
	rfsourcefile = "n/a";
	epss = 0.123;
	epsk = 0.124;
	maxit = 999;
	keff = .333;
	numtypes = 69;
	numel = 381;
	numgroups = 33;
	spaceord = 314;
	angleord = -4;
	faces = 5;
	degfree = 1;
    // newton solver things
    pctype = 1;    // default ilu(0)
    ilulevel = 0; // 
    //
	printout = -1;
	plotout = -1;
    outfile  = "n/a";
	fluxfile = "n/a";
	//
	bcl = -1;
	bcr = -1;
	bcb = -1;
	bct = -1;
    bcn = -1;
    bcs = -1;
	elemx = -2;
	elemy = -2;
	elemz = -2;
}

//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
GlobalInput::~GlobalInput()
{
	return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief This processes the global input.
 *
 * Warning! There is a substantial amount of error checking that should be
 * implemented.  One example is ensuring the matrix of element placements
 * contains the appropriate numbers of columns and rows and zplanes.  There
 * are likely other examples.
 * 
 * Also, it might be worth modifying the xml grammar to specify the underlying
 * data types (e.g. integer, double, etc.).
 *
 */
bool GlobalInput::processFile()
{

    //--------------------------------------------------------------------------
	// get problem name
    xmlNodePtr cur = getXpathElement("/globalinput");
    this->file = (getXpathContent(cur,"name"));
	cout << " inp.name = " << this->file << endl;

    //--------------------------------------------------------------------------
	// move on to control
    cur      = getXpathElement("/*/globalcontrol");
    rfsource = atoi( ( getXpathContent(cur,"rfsource") ) );
    rfsink   = atoi( ( getXpathContent(cur,"rfsink") ) );
    if ( nodeExists( cur, "rfsourcefile" ) )
    	rfsourcefile = getXpathContent( cur, "rfsourcefile" );
    if ( nodeExists( cur, "dbfile" ) )
    	dbfile = getXpathContent( cur, "dbfile" );
    if ( nodeExists( cur, "epsk" ) )	
    	epsk = atof((getXpathContent(cur,"epsk")));
    if ( nodeExists( cur, "epss" ) )	
    	epss = atof((getXpathContent(cur,"epss")));
    if ( nodeExists( cur, "maxit" ) )	
    	maxit = atoi((getXpathContent(cur,"maxit")));
    if ( nodeExists( cur, "keff" ) )	
    	keff = atof((getXpathContent(cur,"keff")));
    numtypes  = atoi((getXpathContent(cur,"numtypes")));
    numel     = atoi((getXpathContent(cur,"numel")));
    numgroups = atoi((getXpathContent(cur,"numgroups")));
    spaceord  = atoi((getXpathContent(cur,"spaceord")));
    angleord  = atoi((getXpathContent(cur,"angleord")));
    faces     = atoi((getXpathContent(cur,"faces")));
	degfree    = numgroups*(1+spaceord)*(1+angleord)*faces*numel;
    if ( nodeExists( cur, "pctype" ) )	
        pctype     = atoi((getXpathContent(cur,"pctype")));
    if ( nodeExists( cur, "ilulevel" ) )	
        ilulevel   = atoi((getXpathContent(cur,"ilulevel")));
    //--------------------------------------------------------------------------
	// move on to output
    cur      = getXpathElement("/*/globaloutput");
    printout = atoi((getXpathContent(cur,"printout")));
    plotout  = atof((getXpathContent(cur,"plotout")));
    if ( nodeExists( cur, "outfile" ) )	
    	outfile = getXpathContent(cur,"outfile");
    if ( nodeExists( cur, "fluxfile" ) )	
    	fluxfile = getXpathContent(cur,"fluxfile");

    //--------------------------------------------------------------------------
	// move on to geometry
    cur      = getXpathElement("/*/globalgeometry");
    bcl = atoi((getXpathContent(cur,"bcl")));
    bcr = atoi((getXpathContent(cur,"bcr")));
    if ( nodeExists( cur, "bcb" ) )	
    	bcb = atoi((getXpathContent(cur,"bcb")));
    if ( nodeExists( cur, "bct" ) )	
    	bct = atoi((getXpathContent(cur,"bct")));
    if ( nodeExists( cur, "bcn" ) )	
    	bcn = atoi((getXpathContent(cur,"bcn")));
    if ( nodeExists( cur, "bcs" ) )	
    	bcs = atoi((getXpathContent(cur,"bcs")));
    elemx = atoi((getXpathContent(cur,"elemx")));
    if ( nodeExists( cur, "elemy" ) )	
    	elemy = atoi((getXpathContent(cur,"elemy")));
    if ( nodeExists( cur, "elemz" ) )	
    	elemz = atoi((getXpathContent(cur,"elemz")));
    
    // try counting the z planes
//    cur = getXpathElement("/*/globalgeometry");
//    getXpathContent(cur,"zplane");
    cur = NULL;

    elements.resize(elemx); // elements[elemx][elemy][elemz]
    for (int i = 0; i < elemx; ++i)
	{
        elements[i].resize(elemy);
        for (int j = 0; j < elemy; ++j)
            elements[i][j].resize(elemz);
	}

    // xmlNodeSetPtr->nodeTab is an array of xmlNodePtrs to the results. 
    // xmlNodesetPtr->nodeNr is the length of the ->nodeTab array.

    xmlNodeSetPtr zplanes = getXpathElements("/*/globalgeometry/zplane");
	xmlNodeSetPtr rows;
	xmlNodeSetPtr cols;
	int i,j; // elements(i,j,k)
	xmlNodePtr cur2;
	xmlChar *key;
	for ( int k = 0; k < zplanes->nodeNr; k++ )
	{
//        cout << " zplane " << k << endl;
//    	cout << "znum = " << getXpathContent(zplanes->nodeTab[k],"number") << endl;
//        cout << " node name = " << zplanes->nodeTab[k]->name << endl;
 	    cur = zplanes->nodeTab[k]->xmlChildrenNode;
		i = 0;
        while ( cur != NULL )
		{
			if ((!xmlStrcmp(cur->name, (const xmlChar *)"r")))
			{
				// now do what
 	    		cur2 = cur->xmlChildrenNode;
				j = 0;
				// loop through columns and get their value
				while ( cur2 != NULL )
				{
//					cout  << " cur2->name = " << cur2->name << endl;
					if ((!xmlStrcmp(cur2->name, (const xmlChar *)"c")))
					{
//						cout << "[i,j,k]=["<<i<<","<<j<<","<<k<<"] "<<(const char*)cur2->xmlChildrenNode->content<<endl;
						elements[i][j][k] = atoi((const char*)cur2->xmlChildrenNode->content)-1;
						j++;
					}		
					cur2 = cur2->next;
				} // while cur2

				i++;
			}
			cur = cur->next;
		} //  while cur
	}

	return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief This echos the input.
 *
 */
void GlobalInput::echoInput()
{
    cout << endl;
    cout << "Echo of user input... " << endl;
    cout << endl;
	cout << " Problem name:  " << file << endl;
    cout << " --------------------------------------------- " << endl;
	cout << " Global Control" << endl;
    cout << " --------------------------------------------- " << endl;
    cout << "       rfsource: " << rfsource << endl;
    cout << "         rfsink: " << rfsink << endl;
    cout << "   rfsourcefile: " << rfsourcefile << endl;
    cout << "         dbfile: " << dbfile << endl;
    cout << "           epss: " << epss << endl;
    cout << "           epsk: " << epsk << endl;
    cout << "          maxit: " << maxit << endl;
    cout << "           keff: " << keff << endl;
    cout << "       numtypes: " << numtypes << endl;
    cout << "          numel: " << numel << endl;
    cout << "      numgroups: " << numgroups << endl;
    cout << "       spaceord: " << spaceord << endl;
    cout << "       angleord: " << angleord << endl;
    cout << "          faces: " << faces << endl;
	cout << "        degfree: " << degfree << endl;
    cout << " --------------------------------------------- " << endl;
	cout << " Global Output" << endl;
    cout << " --------------------------------------------- " << endl;
    cout << "       printout: " << printout << endl;
    cout << "        plotout: " << plotout << endl;
    cout << "        outfile: " << outfile << endl;
    cout << "       fluxfile: " << fluxfile << endl;
    cout << " --------------------------------------------- " << endl;
	cout << " Global Geometry" << endl;
    cout << " --------------------------------------------- " << endl;
    cout << "            bcl: " << bcl << endl;
    cout << "            bcr: " << bcr << endl;
    cout << "            bcb: " << bcb << endl;
    cout << "            bct: " << bct << endl;
    cout << "            bcn: " << bcn << endl;
    cout << "            bcs: " << bcs << endl;
    cout << "          elemx: " << elemx << endl;
    cout << "          elemy: " << elemy << endl;
    cout << "          elemz: " << elemz << endl;
    // verify elements
    for ( int k = 0; k < elemz; ++k )
    {
        cout << "        z plane: " << k << " --- " << endl;
        for(int i = 0; i < elemx; ++i)
        {
			cout << "                 ";
            for(int j = 0; j < elemy; ++j)
		    {
			    cout << elements[i][j][k]+1<< " ";
		    }
            cout << endl;
          
        }
    } // end loop over elements

	return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief This gives a code announcement.
 *
 */
void GlobalInput::sayHello()
{
    time_t t;
    time(&t);
    std::cout << std::ctime(&t);
    cout << " " << endl;
    cout << "|------------------------------------------------|" << endl;
    cout << "|   _____ ______________  ________ _   _ _____   |" << endl;
    cout << "|  /  ___|  ___| ___ \\  \\/  |  ___| \\ | |_   _|  |" << endl;
    cout << "|  \\ `--.| |__ | |_/ / .  . | |__ |  \\| | | |    |" << endl;
    cout << "|   `--. \\  __||    /| |\\/| |  __|| . ` | | |    |" << endl;
    cout << "|  /\\__/ / |___| |\\ \\| |  | | |___| |\\  | | |    |" << endl;
    cout << "|  \\____/\\____/\\_| \\_\\_|  |_|____/\\_| \\_/ \\_/    |" << endl;
    cout << "|                                                |" << endl;
    cout << "|  Solving Eigenvalue Response Matrix Equations  |" << endl;
    cout << "|            using Newton's Technique            |" << endl;
    cout << "|                                                |" << endl;
    cout << "|  Version: 0.1, 04/03/2011                      |" << endl;
    cout << "|------------------------------------------------|" << endl;
    cout << " " << endl;
    cout << " Run on: " << ctime(&t);
}

//---------------------------------------------------------------------------//
//                 end of GlobalInput.cc
//---------------------------------------------------------------------------//

