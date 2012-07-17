//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DBC.hh
 * \author Jeremy Roberts
 * \date   Wed Jan  2 11:31:52 2008
 * \brief  DBC class and type definitions.
 * \note   Copyright (C) Jeremy Roberts, but derived from T. Evan's Denovo.
 */
//---------------------------------------------------------------------------//
// $Rev:: 143                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-09-16 15:40:00 -0400 (Fri, 16 Sep 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef DBC_HH
#define DBC_HH

#include <stdexcept>
#include <string>

// add the serment package configure
#include "serment_config.h"

namespace util
{

//===========================================================================//
/*!
 * \class assertion Exception notification class for Nemesis specific
 * assertions.
 *
 * This class is derived from std::runtime_error.  In fact, this class
 * provides no significant change in functionality from std::runtime_error.
 * This class provides the following features in addition to those found in
 * std::runtime_error: 
 *
 * -# nemesis::assertion does provide an alternate constructor that allows us
 *    to automatically insert file name and line location into error messages.  
 * -# It provides a specialized form of std::runtime_error.  This allows
 *    nemesis code to handle nemesis specific assertions differently from
 *    generic C++ or STL exceptions.  For example
 *    \code
 *    try
 *    {
 *       throw nemesis::assertion( "My error message." );
 *    } 
 *    catch ( nemesis::assertion &a ) 
 *    {
 *       // Catch nemesis exceptions first.
 *       cout << a.what() << endl;
 *       exit(1);
 *    }
 *    catch ( std::runtime_error &e )
 *    { 
 *       // Catch general runtime_errors next
 *       cout << e.what() << endl;
 *    }
 *    catch ( ... )
 *    {
 *       // Catch everything else
 *        exit(1);
 *    }
 *    \endcode
 *
 * \note Assertion should always be thrown as objects on the stack and caught
 *       as references.
 *
 * \sa \ref Nemesis_DBC 
 */
/*!
 * \example harness/test/tstAssert.cc
 * 
 * Assertion and DBC examples.
 */
//===========================================================================//

class assertion : public std::logic_error
{
  public:
    /*!
     * \brief Default constructor for ds++/assertion class.
     *
     * This constructor creates a ds++ exception object.  This object is
     * derived form std::runtime_error and has identical functionality.  The
     * principal purpose for this class is to provide an exception class that
     * is specialized for nemesis.  See the notes on the overall class for more
     * details.
     *
     * \param msg The error message saved with the exception.
     */
    explicit assertion( std::string const & msg )
	:  std::logic_error( msg )
    { /* empty */ }

    /*!
     * \brief Specialized constructor for nemesis::assertion class.
     *
     * This constructor creates a ds++ exception object.  This object is
     * derived form std::runtime_error and has identical functionality.  This
     * constructor is specialized for use by nemesis DbC commands (Require,
     * Ensure, Check, and Insist).  It forms the error message from the test
     * condition and the file and line number of the DbC command.
     *
     * \param cond The expression that failed a DbC test.
     * \param file The source code file name that contains the DbC test.
     * \param line The source code line number that contains the DbC test.
     */
    assertion( std::string const & cond, 
	       std::string const & file, 
	       int const line )
	: std::logic_error( build_message( cond, file, line ) )
    { /* empty */ }

    /*! \brief Destructor for ds++/assertion class.
     * We do not allow the destructor to throw! */
    virtual ~assertion() throw() { /* empty */ }

  private:
    //! Helper function to build error message that includes source file name and line number.
    std::string build_message( std::string const & cond, 
			       std::string const & file, 
			       int         const   line ) const;
};

//---------------------------------------------------------------------------//
// FREE NAMESPACE FUNCTIONS
//---------------------------------------------------------------------------//

//! Throw a nemesis::assertion for Require, Check, Ensure.
void toss_cookies( std::string const & cond, 
		   std::string const & file, 
		   int         const line );

void toss_cookies_ptr(char const * const cond,
		      char const * const file, 
		      int         const line );

//! Throw a nemesis::assertion for Insist.
void insist( std::string const & cond, 
	     std::string const & msg, 
	     std::string const & file, 
	     int         const line);

//! Pointer version of insist
void insist_ptr(char const * const cond, 
		char const * const msg, 
		char const * const file, 
		int          const line);

} // end namespace nemesis

//---------------------------------------------------------------------------//
/*!
 * \page Nemesis_DBC Using the Nemesis Design-by-Contract Macros
 *
 * \section ddbc Using the Nemesis Design-by-Contract Macros
 *
 * The assertion macros are intended to be used for validating preconditions
 * which must be true in order for following code to be correct, etc.  For
 * example, 
 * 
 * \code 
 * Assert( x > 0. ); 
 * y = sqrt(x); 
 * \endcode
 * 
 * If the assertion fails, the code should just bomb.  Philosophically, it
 * should be used to feret out bugs in preceding code, making sure that prior
 * results are within reasonable bounds before proceeding to use those
 * results in further computation, etc.
 * 
 * These macros are provided to support the Design By Contract formalism.
 * The activation of each macro is keyed off a bit in the DBC macro which can 
 * be specified on the command line:
\verbatim
     Bit     DBC macro affected
     ---     ------------------
      0      Require
      1      Check
      2      Ensure
\endverbatim
 * 
 * So for instance, \c -DDBC=7 turns them all on, \c -DDBC=0 turns them all
 * off, and \c -DDBC=1 turns on \c Require but turns off \c Check and \c
 * Ensure.  The default is to have them all enabled.
 *
 * The \c Insist macro is akin to the \c Assert macro, but it provides the
 * opportunity to specify an instructive message.  The idea here is that you
 * should use Insist for checking things which are more or less under user
 * control.  If the user makes a poor choice, we "insist" that it be
 * corrected, providing a corrective hint.
 * 
 * \note We provide a way to eliminate assertions, but not insistings.  The
 * idea is that \c Assert is used to perform sanity checks during program
 * development, which you might want to eliminate during production runs for
 * performance sake.  Insist is used for things which really really must be
 * true, such as "the file must've been opened", etc.  So, use \c Assert for
 * things which you want taken out of production codes (like, the check might
 * inhibit inlining or something like that), but use Insist for those things
 * you want checked even in a production code.
 */
/*!
 * \def Require(condition)
 * 
 * Pre-condition checking macro.  On when DBC & 1 is true.
 */
/*!
 * \def Check(condition)
 * 
 * Intra-scope checking macro.  On when DBC & 2 is true.
 */
/*!
 * \def Ensure(condition)
 * 
 * Post-condition checking macro.  On when DBC & 4 is true.
 */
/*!
 * \def Remember(code)
 * 
 * Add code to compilable code.  Used in the following manner:
 * \code
 *     Remember (int old = x;)
 *     // ...
 *     Ensure (x == old);
 * \endcode
 * On when DBC & 4 is true.
 */
/*!
 * \def Insist(condition, message)
 * 
 * Inviolate check macro.  Insist is always on.
 */
/*!
 * \def Insist_ptr(condition, message)
 * 
 * Same as Insist, except that it uses char pointers, rather than strings.
 * This is more efficient when inlined.
 */
//---------------------------------------------------------------------------//

#ifdef SERMENT_DBC
#define REQUIRE_ON
#define CHECK_ON
#define REMEMBER_ON
#endif

#ifdef REQUIRE_ON
#define Require(c) if (!(c)) util::toss_cookies( #c, __FILE__, __LINE__ )
#else
#define Require(c) 
#endif

#ifdef CHECK_ON
#define Check(c) if (!(c)) util::toss_cookies( #c, __FILE__, __LINE__ )
#define Assert(c) if (!(c)) util::toss_cookies( #c, __FILE__, __LINE__ )
#else
#define Check(c) 
#define Assert(c) 
#endif

#ifdef REMEMBER_ON
#define Ensure(c) if (!(c)) util::toss_cookies( #c, __FILE__, __LINE__ )
#define Remember(c) c
#else
#define Ensure(c) 
#define Remember(c)
#endif

#define Insist(c,m) if (!(c)) util::insist( #c, m, __FILE__, __LINE__ )
#define Insist_ptr(c,m) if (!(c)) util::insist_ptr( #c, m, __FILE__, __LINE__ )

#endif // harness_DBC_hh

//---------------------------------------------------------------------------//
//              end of DBC.hh
//---------------------------------------------------------------------------//
