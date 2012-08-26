//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseIndex.hh
 * \brief  ResponseIndex class definition
 * \author Jeremy Roberts
 * \date   Aug 24, 2012
 */
//---------------------------------------------------------------------------//

#ifndef RESPONSEINDEX_HH_
#define RESPONSEINDEX_HH_

namespace erme_response
{

/*!
 *  \struct ResponseIndex
 *  \brief  Convenience container for indices
 */
struct ResponseIndex
{
  typedef unsigned int size_type;

  ResponseIndex(size_type s  = 0,
                size_type e  = 0,
                size_type p  = 0,
                size_type a  = 0,
                size_type s0 = 0,
                size_type s1 = 0,
                bool      eo = false,
                size_type l = 0)
  : surface(s), energy(e), polar(p), azimuth(a),
    space0(s0), space1(s1), even_odd(eo), local(l)
  {/* ... */}

  size_type surface;
  size_type energy;
  size_type polar;
  size_type azimuth;
  size_type space0;
  size_type space1;
  bool      even_odd;
  size_type local;
};

} // end namespace erme_response

#endif // RESPONSEINDEX_HH_ 

//---------------------------------------------------------------------------//
//              end of file ResponseIndex.hh
//---------------------------------------------------------------------------//
