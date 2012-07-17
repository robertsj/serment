//----------------------------------*-C++-*----------------------------------//
/*
 * \file   InputXML.hh
 * \author Jeremy Roberts
 * \date   04/10/2011
 * \brief  A general abstract input class using XML
 * \note   Inspired by K. Huff's input class from cyclus
 */ 
//---------------------------------------------------------------------------//
// $Rev:: 106                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-06-15 20:35:53 -0400 (Wed, 15 Jun 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef INPUTXML_HH
#define INPUTXML_HH

#include <iostream>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>
#include <libxml/relaxng.h>

using namespace std;

//===========================================================================//
/*
 * \struct xmlFileInfo
 * \brief Hold all information about current XML file
 *
 */ 
//===========================================================================//
struct xmlFileInfo 
{
    /// filename currently being processed
    string filename;  
    /// pointer to the schema URI used for this file
    string* schema;   
    /// XML doc ptr for input file
    xmlDocPtr doc;    
    /// XML XPath search context
    xmlXPathContextPtr xpathCtxt; 
}; 

//===========================================================================//
/*
 * \class InputXML
 * \brief An abstract class with several xml handling methods.
 *
 *  This class uses the libxml2 library to parse xml inputs.  Some of the code
 *  is very much inspired by K. Huff's cyclus input processing routines.
 *
 */ 
//===========================================================================//
class InputXML
{
  	public:

		//----------------------------------------------------------------------		
		// public methods

		InputXML(string S);

		~InputXML();

		bool readInput( char *filename );

	private:

		//----------------------------------------------------------------------
		// private variables

		/// \brief Schema for validating input
		string Schema;

		//----------------------------------------------------------------------		
		// private methods

		void loadFile( char *filename );

		xmlDocPtr validateFile( xmlFileInfo *fileInfo );

	protected:

		//----------------------------------------------------------------------
        // protected variables

		/// \brief Pointer to current XML file
		xmlFileInfo *curFilePtr;

		//----------------------------------------------------------------------		
		// protected methods

		virtual bool processFile() = 0;

		xmlNodePtr getXpathElement( const char* expression );

		xmlNodePtr getXpathElement( xmlNodePtr cur, const char* expression );

		xmlNodeSetPtr getXpathElements( const char* expression );

		xmlNodeSetPtr getXpathElements( xmlNodePtr current, const char* expression );

		const char* getXpathContent( xmlNodePtr cur, const char* expression );

		bool nodeExists( xmlNodePtr current, const char* expression );

        virtual void sayHello()
		{ 
			cout << " YOU SHOULD REDEFINE sayHello! " << endl; 
		};

};

#endif // INPUTXML_HH

//---------------------------------------------------------------------------//
//                 end of InputXML.hh
//---------------------------------------------------------------------------//

