//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DBC.cc
 * \author Thomas M. Evans
 * \date   Wed Jan  2 11:31:52 2008
 * \brief  DBC member definitions.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Rev:: 106                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-06-15 20:35:53 -0400 (Wed, 15 Jun 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#include <sstream>
#include "DBC.hh"

namespace util
{

//===========================================================================//
// ASSERTION CLASS MEMBERS
//===========================================================================//
/*! 
 * Build the error string (private member function).
 * \param cond Condition (test) that failed.
 * \param file The name of the file where the assertion was tested.
 * \param line The line number in the file where the assertion was tested.
 * \retval myMessage A string that contains the failed condition, the file and
 * the line number of the error.
 */
std::string assertion::build_message( std::string const & cond, 
				      std::string const & file, 
				      int         const line ) const
{
    std::ostringstream myMessage;
    myMessage << "Assertion: "
	      << cond
	      << ", failed in "
	      << file
	      << ", line "
	      << line
	      << "." << std::endl;
    return myMessage.str();
}

//===========================================================================//
// FREE FUNCTIONS
//===========================================================================//
/*!
 * \brief Throw a nemesis::assertion for Require, Check, Ensure macros.
 * \return Throws an assertion.
 * \note We do not provide unit tests for functions whose purpose is to throw
 * or exit.
 */
void toss_cookies( std::string const & cond, 
		   std::string const & file, 
		   int         const line )
{
    throw assertion( cond, file, line );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Throw a nemesis::assertion for Require, Check, Ensure macros.
 * \return Throws an assertion.
 * \note We do not provide unit tests for functions whose purpose is to throw
 * or exit.
 */
void 
toss_cookies_ptr(char const * const cond, 
		 char const * const file, 
		 int  const line )
{
    throw assertion( cond, file, line );
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Throw a nemesis::assertion for Insist macros.
 */
void insist( std::string const & cond, 
	     std::string const & msg, 
	     std::string const & file, 
	     int         const line)
{
    std::ostringstream myMessage;
    myMessage <<  "Insist: " << cond << ", failed in "
	      << file << ", line " << line << "." << std::endl
	      << "The following message was provided:" << std::endl
	      << "\"" << msg << "\"" << std::endl;
    throw assertion( myMessage.str() );
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Throw a nemesis::assertion for Insist_ptr macros.
 *
 * Having a (non-inlined) version that takes pointers prevents the compiler
 * from having to construct std::strings from the pointers each time.  This is
 * particularly important for things like SP::operator->, that (a) have an
 * insist in them, (b) don't need complicated strings and (c) are called
 * frequently.
 */
void insist_ptr( char const * const cond, 
		 char const * const msg, 
		 char const * const file, 
		 int          const line)
{
    // Call the other insist for consistency
    insist(cond, msg, file, line);
}

} // end namespace util

//---------------------------------------------------------------------------//
//                 end of DBC.cc
//---------------------------------------------------------------------------//
