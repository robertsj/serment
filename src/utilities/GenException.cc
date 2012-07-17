//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   GenException.cc
 * \author Jeremy Roberts
 * \date   04/09/2011
 * \brief  Member definitions of class GenException
 * \note   Modified version of K. Huff's class from cyclus
 */
//---------------------------------------------------------------------------//
// $Rev:: 106                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-06-15 20:35:53 -0400 (Wed, 15 Jun 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#include <sstream>
#include "GenException.hh"

//namespace util
//{

string itoa(int i)    { stringstream out; out << i; return out.str(); };
string dtoa(double d) { stringstream out; out << d; return out.str(); };

string GenException::prepend = "serment exception";

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
GenException::GenException()
{
	myMessage = prepend;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
GenException::GenException(int line, string file, string msg)
{
	myMessage = prepend + "\n" 
              + "           on line: " + itoa(line) + "\n"
              + "           in file: " + file + "\n" 
              + "           message: " + msg;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
const char* GenException::what() const throw()
{
	//const char* toRet = myMessage;
	//	return toRet;
	return myMessage.c_str();
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
GenException::~GenException() throw() 
{

}

//} // end namespace util
