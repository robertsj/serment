//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   InputXML.cc
 * \author Jeremy Roberts
 * \date   04/10/2011
 * \brief  Member definitions of class InputXML
 * \note   Several functions based on code by K. Huff from the cyclus project.
 */
//---------------------------------------------------------------------------//
// $Rev:: 106                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-06-15 20:35:53 -0400 (Wed, 15 Jun 2011) $:Date of last commit
//---------------------------------------------------------------------------//
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <math.h>
#include "InputXML.hh"
#include "GenException.hh"

using namespace std;

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
InputXML::InputXML(string S)
 : Schema(S)
{
	return;
}

//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
InputXML::~InputXML()
{
	return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief This opens, verifies, and processes an xml input file.
 *
 * The opening and verification are meant to be General enough for use in all
 * input scenarios needing a "load" and "verify" sequence.  The processing is
 * of course program-specific, and will be left as a virtual method in this
 * abstract class.
 *
 * \param filename    the xml input file
 *
 */
bool InputXML::readInput( char *filename )
{

	sayHello();
    cout << " Opening user file: " << filename << endl;

	// create a new xml file pointer
	curFilePtr = new xmlFileInfo;

    // open and validate file, and exit on error
    loadFile(filename);

    // process file
	if ( !processFile() ){
		throw GenException(__LINE__,__FILE__,
                                 "Error processing input.");
		return false;
	}
		
	// delete the xml file point
    delete curFilePtr;

    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Load the XML input file.
 *
 */
void InputXML::loadFile( char *filename )
{

	xmlFileInfo &inputFile = *curFilePtr;

	// Validate the input file agains MainSchema
	inputFile.filename = filename;
  	inputFile.schema = &this->Schema;
  	inputFile.doc = validateFile(&inputFile);

  	// Create xpath evaluation context
  	inputFile.xpathCtxt = xmlXPathNewContext(inputFile.doc);
  	if(inputFile.xpathCtxt == NULL) 
	{
		throw GenException(__LINE__,__FILE__,
                           "Error loading file.");
  	}

	return;
} 

//---------------------------------------------------------------------------//
/*!
 * \brief Validate the XML input.
 *
 */
xmlDocPtr InputXML::validateFile( xmlFileInfo *fileInfo )
{
	xmlRelaxNGParserCtxtPtr ctxt = xmlRelaxNGNewParserCtxt(fileInfo->schema->c_str());
  	if ( NULL == ctxt )
  	{
		throw GenException(__LINE__,__FILE__,
                           "Failed to Generate parser from schema!");
  	}
	xmlRelaxNGPtr schema = xmlRelaxNGParse(ctxt);
	xmlRelaxNGValidCtxtPtr vctxt = xmlRelaxNGNewValidCtxt(schema);
    xmlDocPtr doc;
	doc = xmlReadFile( fileInfo->filename.c_str(), NULL, 0 );
	if ( NULL == doc ) 
	{
		throw GenException(__LINE__,__FILE__,
                           "Failed to parse " + fileInfo->filename );
  	}
  	if ( xmlRelaxNGValidateDoc( vctxt, doc ) ) 
	{
		cout << " File " << fileInfo->filename << " is NOT valid against schema " 
             << *(fileInfo->schema) << endl;
	}
  	else
	{
    	cout << " File " << fileInfo->filename << " is valid against schema " 
    	     << *(fileInfo->schema) << endl;
	}
  	return doc;	
}

//---------------------------------------------------------------------------//
/*!
 * \brief This gets a pointer to the matching node.
 *
 * \param current         pointer to current node
 * \param expression      the name of the element we seek
 */
xmlNodePtr InputXML::getXpathElement( xmlNodePtr current, 
                                      const char* expression )
{
	xmlXPathContextPtr xpathCtxt = curFilePtr->xpathCtxt;
  	xpathCtxt->node = current;
  	// Evaluate xpath expression
  	xmlXPathObjectPtr xpathObj = 
                  xmlXPathEvalExpression((const xmlChar*)expression, xpathCtxt);
  	if(xpathObj == NULL) 
	{
    	fprintf(stderr,"Error: unable to evaluate xpath expression \"%s\"\n", expression);
    	xmlXPathFreeContext(xpathCtxt); 
  	}
  	return xpathObj->nodesetval->nodeTab[0];
}

//---------------------------------------------------------------------------//
/*!
 * \brief This gets a pointer to the matching node.
 *
 * \param expression      the name of the element we seek
 */
xmlNodePtr InputXML::getXpathElement( const char* expression )
{
    return getXpathElement(curFilePtr->doc->children,expression);
}

//---------------------------------------------------------------------------//
/*!
 * \brief This gets a pointer to an array of matching nodes.
 *
 * \param current         pointer to current node
 * \param expression      the name of the element we seek
 */
xmlNodeSetPtr InputXML::getXpathElements( xmlNodePtr current, 
                                          const char* expression )
{

  xmlXPathContextPtr xpathCtxt = curFilePtr->xpathCtxt;
  xpathCtxt->node = current;
  
  /* Evaluate xpath expression */
  xmlXPathObjectPtr xpathObj = 
                  xmlXPathEvalExpression((const xmlChar*)expression, xpathCtxt);
  if(xpathObj == NULL) {
    fprintf(stderr,"Error: unable to evaluate xpath expression \"%s\"\n", expression);
    xmlXPathFreeContext(xpathCtxt); 
  }

  return xpathObj->nodesetval;

}

//---------------------------------------------------------------------------//
/*!
 * \brief This gets a pointer to an array of matching nodes.
 *
 * \param expression      the name of the element we seek
 */
xmlNodeSetPtr InputXML::getXpathElements(const char* expression)
{
    return getXpathElements(curFilePtr->doc->children,expression);
};

//---------------------------------------------------------------------------//
/*!
 * \brief This extracts the content from a given node.
 *
 * \param current         pointer to current node
 * \param expression      the name of the element we seek
 */
const char* InputXML::getXpathContent( xmlNodePtr current, 
                                       const char* expression )
{
  	xmlXPathContextPtr xpathCtxt = curFilePtr->xpathCtxt;
  	xpathCtxt->node = current;
  	// Evaluate xpath expression
  	xmlXPathObjectPtr xpathObj = 
                  xmlXPathEvalExpression((const xmlChar*)expression, xpathCtxt);
	//fprintf(stderr,"---> xpath expression \"%s\"\n", expression);
  	if(xpathObj == NULL) 
	{
		//throw GenException(__LINE__,__FILE__,
        //                 "unable to evaluate xpath expression: " + expression );
    	fprintf(stderr,"Error: unable to evaluate xpath expression \"%s\"\n", expression);
    	xmlXPathFreeContext(xpathCtxt); 
  	}
    if(xmlXPathNodeSetIsEmpty(xpathObj->nodesetval))
    {
		// use for debugging
		cout << " empty set: " << (const xmlChar*)expression << endl;
    	return "";
    }
	else
	{
		// use for debugging
		//cout << xmlXPathNodeSetGetLength(xpathObj->nodesetval) << endl;
	}
  	return (const char*)(xpathObj->nodesetval->nodeTab[0]->children->content);
}

//---------------------------------------------------------------------------//
/*!
 * \brief This checks to see if the (optional, hopefully) node exists.
 *
 * \param current         pointer to current node
 * \param expression      the name of the element we seek
 */
bool InputXML::nodeExists( xmlNodePtr current, const char* expression )
{
  	xmlXPathContextPtr xpathCtxt = curFilePtr->xpathCtxt;
  	xpathCtxt->node              = current;
  	xmlXPathObjectPtr xpathObj   = 
                  xmlXPathEvalExpression((const xmlChar*)expression, xpathCtxt);
    if( xmlXPathNodeSetIsEmpty(xpathObj->nodesetval) )
    {
		// use this for debugging
		cout << "I checked " 
             << expression << ", and it's empty. Using default." << endl;
		return false;
    }
	else
	{
		return true;
	}
}

//---------------------------------------------------------------------------//
//                 end of InputXML.cc
//---------------------------------------------------------------------------//

