//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  matrix_shell_fixture.hh
 *  @brief Shell matrix fixture
 *  @note  Copyright (C) 2012 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

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
    : MatrixShell(m, n, this)
    , d_B(m, n, vec_int(m, 3), vec_int(m, 1))
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

  void multiply(Vector &v_in, Vector &v_out)
  {
    d_B.multiply(v_in, v_out);
  }

  void multiply_transpose(Vector &v_in, Vector &v_out)
  {
    d_B.multiply_transpose(v_in, v_out);
  }

private:

  // Regular matrix.
  Matrix d_B;

};

} // end namespace detran

#endif // MATRIX_SHELL_FIXTURE_HH_ 

//----------------------------------------------------------------------------//
//              end of file matrix_shell_fixture.hh
//----------------------------------------------------------------------------//
