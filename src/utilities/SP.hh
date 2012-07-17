//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/SP.hh
 * \author Thomas M. Evans
 * \date   Thu Jan  3 10:28:56 2008
 * \brief  Smart-Pointer (SP) class definition.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: SP.hh,v 1.1 2008/05/19 17:43:52 TMEvans Stab $
//---------------------------------------------------------------------------//

#ifndef utils_SP_hh
#define utils_SP_hh

#include "DBC.hh"

namespace util
{
 
//===========================================================================//
/*!
 * \struct SPref
 * 
 * \brief Reference holder struct for SP class.
 */
//===========================================================================//

struct SPref 
{
    //! Number of references.
    int refs;

    //! Constructor
    SPref(int r = 1) : refs(r) { }
};

//===========================================================================//
/*!
 * \class SP
 * 
 * \brief Smart pointer implementation that does reference counting.
 *
 * The smart pointer provides a "safe" encapsulation for a standard C++
 * pointer.  Consider: A function new's an object and return the pointer to
 * that object as its return value.  Now it is the caller's responsibility to
 * free the object.  What if the caller passes the pointer to other objects
 * or functions?  What if it is not known which will be deleted first or
 * last?
 *
 * Instead the function can return a "smart pointer".  This SP class uses
 * reference counting to determine the number of current users of a pointer.
 * Each time an SP goes out of scope, the reference count is decremented.
 * When the last user of a pointer is done, the pointer is freed.
 *
 * Note: I am calling this a "smart pointer", not a "safe pointer".  There
 * are clearly ways you can hose this.  In particular, when you bind an SP<T>
 * to a T*, you yield all rights to the T*.  You'd better not squirrel the
 * bare pointer away somewhere and expect to clandestinely use it in other
 * ways or places--death will be sure to follow.  Consequently then, the
 * safest way to use this smart pointer, is to bind it to the contained
 * pointer and then always use the smart pointer.  Immediately returning the
 * smart pointer as a return value, allowing the original bare pointer to go
 * out of scope never to be seen again, is one good example of how to use
 * this.
 * 
 * One good example of bad usage is assigning the same dumb pointer to
 * multiple SPs.  Consider:
 * \code
 *     SP<Foo> f1;
 *     SP<Foo> f2;
 *     // ...
 *     Foo *f = new Foo;
 *     f1 = f;
 *     // ...
 *     f2 = f; // bad, now f1 and f2 assume they "own" f!
 * \endcode
 * Unfortunately, there is no way to check if another SP owns the dumb
 * pointer that you give to a SP.  This is simply something that needs to be
 * watched by the programmer.
 */
/*! 
 * \example utils/test/tstSP.cc
 *
 * denovo::SP (smart pointer) usage example.
 */
//===========================================================================//

template<class T>
class SP 
{
  private: 
    // >>> DATA

    //! Raw pointer held by smart pointer.
    T *p;

    //! Pointer to reference counter.
    SPref *r;

  private:
    // >>> IMPLEMENTATION

    // Free the pointer.
    inline void free();

    //! All derivatives of SP are friends. 
    template<class X> friend class SP;

  public:
    //! Default constructor.
    SP() : p(0), r(new SPref) { Ensure (r); Ensure (r->refs == 1); }

    // Explicit constructor for type T *.
    inline explicit SP(T *p_in);

    // Explicit constructor for type X *.
    template<class X>
    inline explicit SP(X *px_in);

    // Copy constructor for SP<T>.
    inline SP(const SP<T> &sp_in);

    // Copy constructor for SP<X>.
    template<class X>
    inline SP(const SP<X> &spx_in);

    //! Destructor, memory is released when count goes to zero.
    ~SP() { free(); }

    // Assignment operator for type T *.
    inline SP<T>& operator=(T *p_in);

    // Assignment operator for type X *.
    template<class X>
    inline SP<T>& operator=(X *px_in);

    // Assignment operator for type SP<T>.
    inline SP<T>& operator=(const SP<T> sp_in);

    // Assignment operator for type SP<X>.
    template<class X>
    inline SP<T>& operator=(const SP<X> spx_in);

    //! Access operator.
    T* operator->() const { Require(p); return p; }

    //! Dereference operator.
    T& operator*() const { Require(p); return *p; }

    //! Get the base-class pointer; better know what you are doing.
    T* bp() const { return p; }

    //! Boolean conversion operator.
    operator bool() const { return p != 0; }

    //! Operator not.
    bool operator!() const { return p == 0; }

    //! Equality operator for T*.
    bool operator==(const T *p_in) const { return p == p_in; }

    //! Inequality operator for T*.
    bool operator!=(const T *p_in) const { return p != p_in; }

    //! Equality operator for SP<T>.
    bool operator==(const SP<T> &sp_in) const { return p == sp_in.p; }

    //! Inequality operator for SP<T>.
    bool operator!=(const SP<T> &sp_in) const { return p != sp_in.p; }
};

} // end namespace denovo

//---------------------------------------------------------------------------//
// INLINE AND TEMPLATE MEMBERS
//---------------------------------------------------------------------------//

#include "SP.i.hh"

#endif // utils_SP_hh

//---------------------------------------------------------------------------//
//              end of utils/SP.hh
//---------------------------------------------------------------------------//
