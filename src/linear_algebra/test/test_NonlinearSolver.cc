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
using std::cout;
using std::endl;
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

typedef Matrix::SP_matrix SP_matrix;

void fill_A(SP_matrix A, int lb, int ub, int n, double d = 0.0)
{
  for (int i = lb; i < ub; ++i)
  {
    if (i == 0)
    {
      double v[] = { 2.0 + d, -1.0 };
      int r[]    = { 0 };
      int c[]    = { 0, 1 };
      A->insert_values(1, r, 2, c, v);
    }
    else if (i == n)
    {
      double v[] = { -1.0, 2.0 + d };
      int r[]    = { i };
      int c[]    =  { i - 1, i };
      A->insert_values(1, r, 2, c, v);
    }
    else
    {
      double v[] = {-1.0, 2.0 + d, -1.0};
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
    double f_L = 0.5 - 0.5*pow(X.norm(X.L2), 2);
    if (Comm::is_last())
    {
      (*f)[d_m] = f_L;
      L         = (*x)[d_m];
    }
    //COUT("L="<<L);
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

  TestJacobian()
  {
    unsigned int global_count = N, local_start, local_count;
    Comm::partition(global_count, local_start, local_count);
    if (Comm::rank() == Comm::size() - 1) local_count++;
    std::vector<int> nnz(local_count, 3);
    std::vector<int> nnz_od(local_count, 0);
    if (Comm::rank() == Comm::size() - 1)
    {
      nnz[local_count-1]    = local_count;
      nnz_od[local_count-1] = global_count - local_count;
    }
    d_matrix = new Matrix(local_count, local_count, nnz, nnz_od);
  }

  void update(SP_vector x)
  {
    using detran_utilities::range;

    // Size of matrix
    int m = x->local_size();
    if (Comm::rank() ==  Comm::size() - 1) --m;
    Vector X(x->local_size(), 0);
    X.copy(*x);

    double lambda = 0.0;
    if (Comm::is_last()) lambda = (*x)[m];
    Comm::broadcast(&lambda, 1, Comm::size() - 1);

    X.scale(-1); // since we insert -X in J

    // Fill all but the last row and column
    int lb = X.lower_bound();
    int ub = X.upper_bound();
    fill_A(d_matrix, lb, ub, d_matrix->number_global_rows()-1, -lambda);

    // Last column
    {
      Vector::vec_int r   = range<int>(lb, ub);
      int             c[] = {d_matrix->number_global_columns() - 1};
      d_matrix->insert_values(X.local_size(), &r[0], 1, c, &X[0]);
    }

    // Last row
    Vector::SP_vector x_seq = X.collect_on_root(Comm::last());
    if (Comm::is_last())
    {
      Vector::vec_int c   = range<int>(0, d_matrix->number_global_columns());
      int             r[] = {d_matrix->number_global_columns() - 1};
      d_matrix->insert_values(1, r, x_seq->local_size(), &c[0], &(*x_seq)[0]);
    }

    d_matrix->assemble();
  }

};

int test_NonlinearSolver(int argc, char *argv[])
{
  Matrix::SP_matrix A;
  unsigned int M = N, s, m;
  Comm::partition(M, s, m);
  std::vector<int> nnz(m, 3), nnz_od(m, 0);
  A = new Matrix(m, m, nnz, nnz_od);
  fill_A(A, A->lower_bound(), A->upper_bound(), A->number_global_rows()-1);
  A->assemble();

  TestResidual::SP_residual f(new TestResidual(A));
  TestJacobian::SP_jacobian J(new TestJacobian());

  Vector::SP_vector x, z;
  x = new Vector(J->matrix()->number_local_rows(), 1.23);
  z = new Vector(J->matrix()->number_local_rows(), 0.0);
  Vector X0(m, 0.0), X1(*x, m);
  for (int i = 0; i < m; ++i)
    X0[i] = A->lower_bound() + i;
  X0.scale(1.0 / X0.norm(X1.L2));

  // power iterations
  double L = 0.0;
  for (int i = 0; i < 1000; ++i)
  {
    A->multiply(X0, X1); L = X1.norm(X1.L2); X1.scale(1.0 / L);
    A->multiply(X1, X0); L = X0.norm(X0.L2); X0.scale(1.0 / L);
  }
  A->multiply(X0, X1); L = X1.norm(X1.L2); X1.scale(1.0 / L);
  if (Comm::is_last())
  {
    (*x)[m] = L;
  }
  // A->display(x->BINARY, "A.out");


  J->update(x);
  //A->display();
  J->matrix()->display(x->BINARY, "J.out");
  NonlinearSolver solver;
  NonlinearSolver::SP_db db = detran_utilities::InputDB::Create();
  solver.setup(db, f, J, J);
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
