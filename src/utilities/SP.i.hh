//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/SP.i.hh
 * \author Thomas M. Evans
 * \date   Thu Jan  3 10:28:56 2008
 * \brief  Member definitions of class SP.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: SP.i.hh,v 1.1 2008/05/19 17:43:34 TMEvans Stab $
//---------------------------------------------------------------------------//

#ifndef utils_SP_i_hh
#define utils_SP_i_hh

namespace util
{

//---------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//---------------------------------------------------------------------------//
/*!
 * \brief Do equality check with a free pointer.
 */
template<class T>
bool operator==(const T *pt, const SP<T> &sp)
{
    return sp == pt;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Do inequality check with a free pointer.
 */
template<class T>
bool operator!=(const T *pt, const SP<T> &sp)
{
    return sp != pt;
}

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Explicit constructor for type T *.
 *
 * This constructor is used to initialize a SP with a pointer, ie.
 * \code
 *     Foo *f = new Foo;
 *     SP<Foo> spf(f);   // f now owned by spf
 *     SP<Foo> spf2 = f; // error! does not do implicit conversion
 * \endcode
 * Once a pointer is "given" to a SP, the SP takes control. This means that
 * spf will delete f when the last SP to f is destroyed.
 *
 * \param p_in pointer to type T
 */
template<class T>
SP<T>::SP(T *p_in)
    : p(p_in),
      r(new SPref)
{
    Ensure (r);
    Ensure (r->refs == 1);
}

//---------------------------------------------------------------------------//
/*!
 * \brief  Explicit constructor for type X *.
 *
 * This constructor is used to initialize a base class smart pointer of type
 * T with a derived class pointer of type X or, equivalently, any types in
 * which X * is convertible to T * through a dynamic_cast.  Consider,
 * \code
 *     class Base {//...}; 
 *     class Derived : public Base {//...};
 *     
 *     SP<Base> spb(new Derived); // spb base class SP to Derived type
 *     Derived *d = new Derived;
 *     SP<Base> spb2(d);          // different syntax
 *     SP<Base> spb3 = d;         // error! no implicit conversion
 *     
 *     Derived *d2;
 *     SP<Base> spb4(d2);         // error! cannot initialize with NULL
 *                                // pointer of different type than T
 * \endcode
 * The pointer to X must not be equal to NULL.  The SP owns the pointer when
 * it is constructed.
 *
 * \param px_in pointer to type X that is convertible to T *
 */
template<class T>
template<class X>
SP<T>::SP(X *px_in)
{
    Require (px_in);

    // make a dynamic cast to check that we can cast between X and T
    T *np = dynamic_cast<T *>(px_in);

    // check that we have made a successfull cast if px exists
    Insist(np, "Incompatible dumb pointer conversion between X and SP<T>.");

    // assign the pointer and reference
    p = np;
    r = new SPref;

    Ensure (p);
    Ensure (r);
    Ensure (r->refs == 1);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Copy constructor for SP<T>.
 *
 * \param sp_in smart pointer of type SP<T>
 */
template<class T>
SP<T>::SP(const SP<T> &sp_in)
    : p(sp_in.p),
      r(sp_in.r)
{
    Require (r);

    // advance the reference to T
    r->refs++;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Copy constructor for SP<X>.
 *
 * This copy constructor requires that X * is convertible to T * through a
 * dynamic_cast.  The pointer in spx_in can point to NULL; however, it must
 * still be convertible to T * through a dynamic_cast.
 *
 * \param spx_in smart pointer of type SP<X>
 */
template<class T>
template<class X>
SP<T>::SP(const SP<X> &spx_in)
{
    Require (spx_in.r);

    // make a pointer to T *
    T *np = dynamic_cast<T *>(spx_in.p);
    Insist(spx_in.p ? np != 0 : true, 
	   "Incompatible SP conversion between SP<X> and SP<T>.");

    // assign the pointer and reference
    p = np;
    r = spx_in.r;
    
    // advance the reference to T
    r->refs++;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment operator for type T *.
 *
 * The assignment operator checks the existing count of the smart pointer and
 * then assigns the smart pointer to the raw pointer.  As in copy
 * construction, the smart pointer owns the pointer after assignment. Here is
 * an example of usage:
 * \code
 *     SP<Foo> f;         // has reference to NULL pointer
 *     f      = new Foo;  // now has 1 count of Foo
 *     Foo *g = new Foo; 
 *     f      = g;        // f's original pointer to Foo is deleted
 *                        // because count goes to zero; f now has
 *                        // 1 reference to g
 * \endcode
 * 
 * \param p_in pointer to T
 */
template<class T>
SP<T>& SP<T>::operator=(T *p_in)
{
    // check if we already own this pointer
    if (p == p_in)
	return *this;

    // next free the existing pointer
    free();
    
    // now make add p_in to this pointer and make a new reference to it
    p = p_in;
    r = new SPref;

    Ensure (r->refs == 1);
    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment operator for type X *.
 *
 * This assignment requires that X * is convertible to T * through a
 * dynamic_cast.  It follows the same principle as SP(X*); however,
 * this is assignment:
 * \code
 *     SP<Base> b;
 *     b = new Derived;
 * \endcode
 * The pointer to X must not be equal to NULL.
 *
 * \param px_in pointer to type X * that is convertible to type T * through a
 * dynamic cast
 */
template<class T>
template<class X>
SP<T>& SP<T>::operator=(X *px_in)
{
    Require (px_in);

    // do a dynamic cast to ensure convertiblility between T* and X*
    T *np = dynamic_cast<T *>(px_in);
    Insist(np, "Incompatible dumb pointer conversion between X and SP<T>.");

    // now assign this to np (using previously defined assignment operator)
    *this = np;
    
    Ensure (r->refs == 1);
    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment operator for type SP<T>.
 *
 * \param sp_in smart pointer of type SP<T>
 */
template<class T>
SP<T>& SP<T>::operator=(const SP<T> sp_in)
{
    Require (sp_in.r);

    // see if they are equal
    if (this == &sp_in || p == sp_in.p)
	return *this;

    // free the existing pointer
    free();

    // assign p and r to sp_in
    p = sp_in.p;
    r = sp_in.r;

    // add the reference count and return
    r->refs++;
    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment operator for type SP<X>.
 *
 * This assignment requires that X * is convertible to T * through a
 * dynamic_cast.  The pointer in spx_in can point to NULL; however, it must
 * still be convertible to T * through a dynamic_cast.
 *
 * \param spx_in smart pointer of type SP<X>
 */
template<class T>
template<class X>
SP<T>& SP<T>::operator=(const SP<X> spx_in)
{
    Require (spx_in.r);

    // make a pointer to T *
    T *np = dynamic_cast<T *>(spx_in.p);
    Insist(spx_in.p ? np != 0 : true, 
	   "Incompatible SP conversion between SP<X> and SP<T>.");

    // check to see if we are holding the same pointer (and np is not NULL);
    // to NULL pointers to the same type are defined to be equal by the
    // standard 
    if (p == np && p)
    {
	// if the pointers are the same the reference count better be the
	// same 
	Check (r == spx_in.r);
	return *this;
    }

    // we don't need to worry about the case where p == np and np == NULL
    // because this combination is impossible; if np is NULL then it belongs
    // to a different smart pointer; in other words, if p == np and np ==
    // NULL then r != spx_in.r

    // free the existing pointer and reference
    free();

    // assign new values
    p = np;
    r = spx_in.r;

    // advance the counter and return
    r->refs++;
    return *this;
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION
//---------------------------------------------------------------------------//
/*!
 * \brief Decrement the count and free the pointer if count is zero.
 *
 * Note that it is perfectly acceptable to call delete on a NULL pointer.
 */
template<class T>
void SP<T>::free()
{
    Require (r);
    
    // if the count goes to zero then we free the data
    if (--r->refs == 0)
    {
	delete p;
	delete r;
    }
}

}

#endif // utils_SP_i_hh

//---------------------------------------------------------------------------//
//              end of utils/SP.i.hh
//---------------------------------------------------------------------------//
