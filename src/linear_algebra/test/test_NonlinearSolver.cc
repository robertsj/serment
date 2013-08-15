//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   test_NonlinearSolver.cc
 *  @brief  Test of Matrix class.
 *  @note   Copyright (C) 2013 Jeremy Roberts.
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                  \
        FUNC(test_NonlinearSolver)

#include "utilities/TestDriver.hh"
#include "linear_algebra/NonlinearSolver.hh"
#include "linear_algebra/Matrix.hh"
#include "linear_algebra/LinearAlgebraSetup.hh"
#include "utilities/MathUtilities.hh"
#include <iostream>

using namespace linear_algebra;
using namespace serment_comm;
using namespace detran_test;
using namespace detran_utilities;

using std::cout;
using std::endl;
using detran_utilities::range;

#define COUT(c) cout << c << endl;

int main(int argc, char *argv[])
{
  linear_algebra::initialize(argc, argv, true);
  RUN(argc, argv);
  linear_algebra::finalize();
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

/**
 *  The test problem is using Newton to solve the eigenvalue problem
 *    A*x = L*x
 *  where A is the matrix [2 -1 0 .. ; -1 2 -1 0 0; ... -1 2]
 *  with a dimension of 100.  We want the largest eigenvalue, which is
 *     3.99903256458398
 *  To solve Ax=Lx via Newton, we define the residual as
 *     f_x = Ax-Lx
 *     f_L = 0.5(1-x'x)
 *  leading to the Jacobian
 *     J = | (A-L*I)   -x |
 *         | -x'        0 |
 */

double lambda_0 = 3.99903256458398;
double lambda_1 = 3.99613119426719;
int N = 100;
double ref_10[] =
    { 0.120131165878581, -0.230530019145232,  0.322252701275550,
     -0.387868386059132,  0.422061280946316, -0.422061280946317,
      0.387868386059134, -0.322252701275552,  0.230530019145233,
     -0.120131165878581,  3.91898594722899};

typedef SP<Matrix> SP_matrix;

void fill_A(SP_matrix A)
{
  int lb = A->lower_bound();
  int ub = A->upper_bound();
  for (int i = lb; i < ub; ++i)
  {
    if (i == 0)
    {
      double v[] = { 2.0, -1.0 };
      int r[]    = { 0 };
      int c[]    = { 0, 1 };
      A->insert_values(1, r, 2, c, v);
    }
    else if (i == ub-1)
    {
      double v[] = { -1.0, 2.0 };
      int r[]    = { i };
      int c[]    = { i - 1, i };
      A->insert_values(1, r, 2, c, v);
    }
    else
    {
      double v[] = {-1.0, 2.0, -1.0};
      int r[]    = {i};
      int c[]    = {i-1, i, i+1};
      A->insert_values(1, r, 3, c, v);
    }
  }
}

class TestResidual: public NonlinearResidualBase
{
public:
  TestResidual(SP_matrix A) : d_A(A)
  {
    d_m = d_A->number_local_rows();
  }
  void evaluate(Vector *x, Vector *f)
  {
    Vector X(*x, d_m), F(*f, d_m);
    d_A->multiply(X, F);
    double L   = 0.0;
    double f_L = 0.5 - 0.5 * pow(X.norm(X.L2), 2);
    if (Comm::is_last())
    {
      (*f)[d_m] = f_L;
      L         = (*x)[d_m];
    }
    Comm::broadcast(&L, 1, Comm::last());
    F.add_a_times_x(-L, X);
  }
  Matrix::SP_matrix A() {return d_A;}
private:
  Matrix::SP_matrix d_A;
  int d_m;
};


class TestJacobian: public JacobianBase
{
public:

  TestJacobian(SP_matrix A, bool as_pc)
    : d_A(A)
    , d_m(d_A->number_local_rows())
    , d_local_size(d_A->number_local_rows())
    , d_as_pc(as_pc)
  {
    if (Comm::is_last())
      ++d_local_size;
    std::vector<int> nnz(d_local_size, 3);
    std::vector<int> nnz_od(d_local_size, 0);
    if (Comm::is_last())
    {
      nnz[d_local_size-1]    = d_local_size;
      nnz_od[d_local_size-1] = d_A->number_global_rows() + 1 - d_local_size;
    }
    d_matrix = new Matrix(d_local_size, d_local_size, nnz, nnz_od);
  }

  void update(SP_vector x)
  {
    Matrix &J = *(dynamic_cast<Matrix*>(d_matrix.bp()));

    Vector X(x->local_size(), 0);
    X.copy(*x);

    double lambda = 0.0;
    if (Comm::is_last())
      lambda = (*x)[d_m];
    Comm::broadcast(&lambda, 1, Comm::last());

    X.scale(-1);

    // Fill all but the last row and column
    J.insert_values(d_A, J.INSERT);
    J.assemble(J.FLUSH);
    for (int i = d_A->lower_bound(); i < d_A->upper_bound(); ++i)
    {
      double val = -lambda;
      J.insert_values(1, &i, 1, &i, &val, J.ADD);
    }
    J.assemble(J.FLUSH);

    // Last column
    {
      Vector::vec_int r   = range<int>(d_A->lower_bound(), d_A->upper_bound());
      int             c[] = {J.number_global_columns() - 1};
      J.insert_values(d_m, &r[0], 1, c, &X[0], J.INSERT);
    }

    // Last row
    Vector::SP_vector x_1 = X.collect_on_root(Comm::last());
    if (Comm::is_last())
    {
      Vector::vec_int c   = range<int>(0, J.number_global_columns());
      int             r[] = {J.number_global_columns() - 1};
      if (d_as_pc)
        (*x_1)[d_m] = 1.0;
      else
        (*x_1)[d_m] = 0.0;
      J.insert_values(1, r, x_1->local_size(), &c[0], &(*x_1)[0], J.INSERT);
    }

    J.assemble();
  }
private:
  SP<Matrix> d_A;
  int d_m;                // local size of A
  int d_local_size;       // # unknowns
  bool d_as_pc;           // flag indicating whether or not to place 1 in corner
};

int test_NonlinearSolver(int argc, char *argv[])
{
  Matrix::SP_matrix A;
  unsigned int M = N, s, m;
  Comm::partition(M, s, m);
  std::vector<int> nnz(m, 3), nnz_od(m, 0);
  A = new Matrix(m, m, nnz, nnz_od);
  fill_A(A);
  A->assemble();


  TestResidual::SP_residual f(new TestResidual(A));
  TestJacobian::SP_jacobian J(new TestJacobian(A, false));
  TestJacobian::SP_jacobian P(new TestJacobian(A, true));

  Vector::SP_vector x, z;
  x = new Vector(J->matrix()->number_local_rows(), 1.23);
  z = new Vector(J->matrix()->number_local_rows(), 0.0);
  Vector X0(m, 0.0), X1(*x, m);
  for (int i = 0; i < m; ++i)
    X0[i] = A->lower_bound() + i;
  X0.scale(1.0 / X0.norm(X1.L2));

  // power iterations
  double L = 0.0;
  for (int i = 0; i < 3000; ++i)
  {
    A->multiply(X0, X1); L = X1.norm(X1.L2); X1.scale(1.0 / L);
    A->multiply(X1, X0); L = X0.norm(X0.L2); X0.scale(1.0 / L);
  }
  A->multiply(X0, X1); L = X1.norm(X1.L2); X1.scale(1.0 / L);
  if (Comm::is_last())
  {
    (*x)[m] = L;
  }
//  x->display(x->BINARY, "X.out");
//  A->display(x->BINARY, "A.out");
//  J->update(x);
//  A->display(x->BINARY, "A.out");
//  J->matrix()->display(x->BINARY, "J.out");

  NonlinearSolver solver;
  NonlinearSolver::SP_db db = detran_utilities::InputDB::Create();
  solver.setup(db, f, J, P);
  solver.solve(x);

  if (Comm::is_last())
  {
    L = (*x)[m];
    COUT("rank = " << Comm::rank() << " eigenvalue " << L)
  }
  Comm::broadcast(&L, 1, Comm::last());
  COUT("L = " << 3.99903256458398 << " " << L)
  TEST(detran_utilities::soft_equiv(L, 3.99903256458398, 1e-9));

  Comm::global_barrier();
  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_NonlinearSolver.cc
//----------------------------------------------------------------------------//
