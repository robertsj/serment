//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mpi_matrix_test.cc
 * \author Jeremy Roberts
 * \date   12/14/2010
 * \brief  Tests building a matrix and vector in parallel, and solves a system.
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 75                                            $:Rev of last commit
// $Author:: bert                                       $:Author of last commit
// $Date:: 2010-12-19 09:47:01 -0500 (Sun, 19 Dec 2010) $:Date of last commit
//---------------------------------------------------------------------------//

#include <iostream>
#include <fstream>
#include "typedefs.hh"
#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"

using namespace std;

//---------------------------------------------------------------------------//
/*!
 * \brief This builds the 1-D finite difference matrix.
 *
 */
void buildmatrix( Mat &A, integer N )
{
    cout << " buildmatrix " << endl;
    MatCreateMPIAIJ( PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                      N-1, N-1, 3, PETSC_NULL, 2, PETSC_NULL, &A);
    integer Istart, Iend;
    MatGetOwnershipRange(A,&Istart,&Iend);
    scalar val[3], val2[2];
    integer col[3], col2[2];
    val[0] = -1.0; val[1] = 2.0; val[2] = -1.0;
    for (integer Ii=Istart; Ii<Iend; Ii++) 
    {
        if ( Ii > 0 and Ii < N-2 )  
        {
            col[0] = Ii-1; col[1] = Ii; col[2] = Ii+1;  
            MatSetValues( A, 1, &Ii, 3, col, val, INSERT_VALUES );
        }
        else if ( Ii == 0 )
        {
            col2[0] = Ii; col2[1] = Ii+1;
            val2[0] = 2.0; val2[1] = -1.0;
            MatSetValues( A, 1, &Ii, 2, col2, val2, INSERT_VALUES );
        }
        else // ( Ii == N-1 )
        {
            col2[0] = Ii-1; col2[1] = Ii;
            val2[0] = -1.0; val2[1] = 2.0;
            MatSetValues( A, 1, &Ii, 2, col2, val2, INSERT_VALUES );
        }
    }
    MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY );
 
    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief This builds the right hand side.
 *
 */
void buildvector( Vec &x, Vec &b, integer N )
{
    VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,N-1,&x);
    VecDuplicate(x,&b);
    scalar h2 = 1.0/N/N;
    integer Istart, Iend;
    VecGetOwnershipRange(b,&Istart,&Iend);
    for (integer Ii=Istart; Ii<Iend; Ii++) 
        VecSetValue( b, Ii, h2, INSERT_VALUES );
    VecAssemblyBegin( b );
    VecAssemblyEnd( b );
    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief This is a simple parallel parallel linear solver demonstration.
 *
 * The system solved is y''(x) = 1 for x = [0,1], subject to y(0)=y(1)=0.  It
 * is implemented as the linear system Ax = b, where A is the 1-D 
 * second-difference operator defined as (1/h)^2 * [ 2 -1 0 ..., -1 2 -1 ....]
 * where h = 1/N, and N is the number of divisions.  The right hand side is
 * a vector of 1's.
 *
 * Uses as a reference: ksp/examples/tutorials/ex2.c
 */
int main(int argc, char *args[])
{
    // initialize petsc
    PetscInitialize(&argc,&args,PETSC_NULL,PETSC_NULL);


    Mat K;      // 1-D finite difference matrix
    Vec x;      // unknowns
    Vec b;      // right hand sid
    integer N=1000000;  // number of divisions: | 1 | 2 | ... | N | --> N-1 unknowns

    integer size;
    integer rank;
    MPI_Comm_size( PETSC_COMM_WORLD, &size );
    MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
    integer print=0;

//    if ( rank==0 )
//    {
//        cout << " input N " << endl;
//        cin >> N;
//        cout << " view things? " << endl;
//        cin >> print;
//    }    
    
    MPI_Bcast( &N, 1, MPI_INT, 0, PETSC_COMM_WORLD );
    MPI_Bcast( &print, 1, MPI_INT, 0, PETSC_COMM_WORLD );

    scalar time = MPI_Wtime();

    buildmatrix( K, N );
    buildvector( x, b, N );

    KSP ksp;
    KSPCreate( PETSC_COMM_SELF, &ksp );
    KSPSetType( ksp, KSPGMRES );
    KSPSetOperators( ksp, K, K, SAME_NONZERO_PATTERN );
    KSPSetFromOptions( ksp );
    KSPSolve( ksp, b, x );

    if (rank==0)
    {
        ofstream myfile;
        myfile.open ("time.txt");
        myfile << " time = " << MPI_Wtime() - time << endl << endl;
        myfile.close();
    }

    if (print==1) VecView( x, PETSC_VIEWER_STDOUT_WORLD );

    MatDestroy( K );
    VecDestroy( b );
    VecDestroy( x );
    KSPDestroy( ksp );

	PetscFinalize();
	cout << "  ... all finished. now check for errors! " << endl;
    return 0;
}


//---------------------------------------------------------------------------//
//                 end of mpi_matrix_test.cc
//---------------------------------------------------------------------------//

