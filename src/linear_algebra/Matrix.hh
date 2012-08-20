//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Matrix.hh
 * \brief  Matrix 
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 */
//---------------------------------------------------------------------------//

#ifndef MATRIX_HH_
#define MATRIX_HH_

#include "MatrixBase.hh"

namespace linear_algebra
{

class Matrix: public MatrixBase
{

public:

  Matrix(const size_type m,
         const size_type n,
         const vec_int &number_nonzeros,
         const vec_int &number_nonzeros_offdiagonal = vec_int(0));

private:

};

} // end namespace linear_algebra

#endif // MATRIX_HH_ 

//---------------------------------------------------------------------------//
//              end of file Matrix.hh
//---------------------------------------------------------------------------//
