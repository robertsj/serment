//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   matrix_shell_fixture.hh
 * \brief  matrix_shell_fixture 
 * \author Jeremy Roberts
 * \date   Aug 21, 2012
 */
//---------------------------------------------------------------------------//

#ifndef MATRIX_SHELL_FIXTURE_HH_
#define MATRIX_SHELL_FIXTURE_HH_

#include "MatrixShell.hh"
#include "Matrix.hh"
#include <iostream>

namespace linear_algebra
{

class MyMatrixShell: public MatrixShell
{

public:

  MyMatrixShell(const size_type m, const size_type n)
  :  MatrixShell(m, n, this)
  ,  d_B(m, n, vec_int(m, 3), vec_int(m, 0))
  {

    for (int row = d_B.lower_bound(); row < d_B.upper_bound(); row++)
    {
      if (row == 0)
      {
        int c[] = {0, 1};
        double v[] = {-2, 1};
        d_B.insert_values(1, &row, 2, c, v);
      }
      else if (row == d_B.number_global_rows() - 1)
      {
        int c[] = {d_B.number_global_rows() - 2, d_B.number_global_rows() - 1};
        double v[] = {1, -2};
        d_B.insert_values(1, &row, 2, c, v);
      }
      else
      {
        int c[] = {row - 1, row, row + 1};
        double v[] = {1, -2, 1};
        d_B.insert_values(1, &row, 3, c, v);
      }
    }
    d_B.assemble();

  }

  PetscErrorCode shell_multiply(Vec x, Vec y)
  {
    return MatMult(d_B.A(), x, y);
  }

  PetscErrorCode shell_multiply_transpose(Vec x, Vec y)
  {
    return MatMultTranspose(d_B.A(), x, y);
  }

private:

  // Regular matrix.
  Matrix d_B;

};

} // end namespace detran

#endif // MATRIX_SHELL_FIXTURE_HH_ 

//---------------------------------------------------------------------------//
//              end of file matrix_shell_fixture.hh
//---------------------------------------------------------------------------//
