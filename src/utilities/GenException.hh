//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   GenException.cc
 * \author Jeremy Roberts
 * \date   04/09/2011
 * \brief  A simple exception class.
 * \note   Modified version of K. Huff's class from cyclus
 */
//---------------------------------------------------------------------------//
// $Rev:: 106                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-06-15 20:35:53 -0400 (Wed, 15 Jun 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef GENEXCEPTION_HH
#define GENEXCEPTION_HH

#include <iostream>
#include <exception>
#include <string>

using namespace std;

//namespace util
//{

/**
 *  A generic mechanism to manually manage exceptions
 */
class GenException: public std::exception
{

protected:

    /// The message associated with this exception.
    std::string myMessage;
    
    /// A string to prepend to all message of this class.
    static std::string prepend;
    
public:
    
    /// Constructs a new GenException with the default message.
    GenException();
    
    /**
     * @brief Constructs a new GenException with a provided message
     *
     * @param line line of code erring
     * @param file file in which error occurs
     * @param msg the message
     */
    GenException(int line, string file, string msg);
    
    /**
     * Returns the error message associated with this GenException.
     *
     * @return the message
     */
    virtual const char* what() const throw();
    
    /**
     * Destroys this GenException.
     */
    virtual ~GenException() throw();
    
};

//} // end namespace util

#endif // GENEXCEPTION_HH
