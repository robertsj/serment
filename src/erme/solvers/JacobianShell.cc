//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   JacobianShell.cc
 * \author Jeremy Roberts
 * \date   10/25/2010
 * \brief  Member definitions of abstract class JacobianShell
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 177                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-12-09 16:31:49 -0500 (Fri, 09 Dec 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#include <iostream>
#include "LinAlg.hh"
#include "GlobalProblem.hh"
#include "Newton.hh"
#include "MyMatVecWrap.hh"
#include "JacobianShell.hh"

//---------------------------------------------------------------------------//
/*!
 * \brief This performs the action of the Jacobian (\f$ \mathbf{f}' \f$ )on a 
 *        vector \f$ \mathbf{f} \f$
 *
 * The resulting vector is \f$ \mathbf{f}'\mathbf{f} \f$, or "fpf" where "p"
 * denotes the prime.
 *
 * Note, the vector \f$ \vec{f} \f$  does not necessarily begin as the
 * residual but rather as a scaled residual. The PETSc manual does not seem
 * to state this explicitly, but a little trial and error confirmed this---so
 * don't fret when debugging!
 *
 * The Jacobian is defined as 
 *  \f[
 *  \mathbf{f'(x)} = \left [\begin{array}{ccc}
 *	  (\mathbf{M}\mathbf{R}-\lambda \mathbf{I})                       &
 *	            \mathbf{M}\mathbf{R_k}\mathbf{J_-}                     &
 *	                      \mathbf{J_-}                                  \\
 *	  (\mathbf{F}-k\mathbf{L})                                        &
 *	            (\mathbf{F_k}-k\mathbf{L_k}-\mathbf{L}) \mathbf{J_-}   &
 *	                      0                                             \\
 *	  \mathbf{J^T_-}                                                  &
 *	            0                                                      &
 *	                      0
 *	\end{array} 
 *  \right ]  \, .
 *  \label{eq:jacobian}
 *  \f]
 * For \f$ \mathbf{R}(k) \f$ of size \f$ m\times m \f$, the Jacobian is 
 * of size \f$ (m+2)\times(m+2) \f$.
 *
 * Most of the Jacobian is known <EM> a priori</EM>, and so most of the action
 * can be computed directly.  Only the first \f$ m-1 \f$ rows of the 
 * \f$ (m-1) \f$th column require approximations via finite differences.
 *
 */
void JacobianShell::myMatVec( Mat &M, Vec &f, Vec &fpf )
{

    bool debug = false;
    Mat JJ;
    if (debug)
    {
        PetscViewer bview;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD,"jac.bin",FILE_MODE_READ ,&bview);
        //MatLoad(bview, MATSEQAIJ, &JJ );
    }

//    cout << " ********** JACOBIAN ********** " << endl;
//    cout << " x is !! " << endl;    
//    VecView( x->V, PETSC_VIEWER_STDOUT_SELF );
//    cout << " f is !! " << endl;    
//    VecView( f, PETSC_VIEWER_STDOUT_SELF );
//    cout << " x is " << endl;
//    x->viewMe();
//    cout << " fpf is !! " << endl;    
//    VecView( fpf, PETSC_VIEWER_STDOUT_SELF );
//    cout << " ******************************* " << endl;

    scalar *x_a, *f_a, *fpf_a; // unknowns x, f(x), and f'(x)*f(x)

    // we need to extract x, f, and fpf
    // recall, x is a SermentVector, so we extract its V
    VecGetArray( x->V, &x_a   );
    VecGetArray( f,    &f_a   );
    VecGetArray( fpf,  &fpf_a );

    // The Jacobian has the form
    //   | (M*R-lambda*I)   M*R_k*Jinc           -Jinc |   
    //   | (F - k*L)        (F_k-k*L_k-L)*Jinc    0    |
    //   | -Jinc'           0                     0    |
    // and so we compute the action in pieces.  The first m-2 entries of the
    // output of f'*f are the sum of the (m-2)x(1) vectors
    //    fpfJ      = M*R*fJ - lambda*fJ +  M*R_k*Jinc*fk + Jinc*flambda
    // the (m-1)th entry is the sum
    //    fpfk      = (F-k*L)*fJ  +  (F_k-k*L_k-L)*Jinc*fk
    // and the final element is
    //    fpflambda = Jinc'*fJ
    // We approximate M*R_k*Jinc as
    //     d/dk(  M*R(k)*J ) ~ M * ( R(2)-R(1) )/eps * J 
    //                       = M*R(2)*J/eps - M*R(1)*J/eps
    // Similarly,
    //     d/dk(F*fJ) = F(2)fJ/eps - F(1)fJ/eps
    //     d/dk(L*fJ) = (A(2)*fJ+L(2)*fJ)/eps - (A(1)*fJ+L(1)*fJ)/eps
    // Note the ambiguous labeling: "L" as written in the Jacobian refers to
    // a total "loss" operator, while in the code "A+L" gives the total 
    // absorption plus leakage, i.e. the total loss.

    // the incident current vector J- comprises the 1st m-2 entries of x
//    scalar  zero      = 0.0;
//    scalar  one       = 1.0;
    scalar  none      = -1.0;
    scalar  fk        = f_a[m-2];  // (m-1)th entry of f
    scalar  flambda   = f_a[m-1];  // (m)th entry of f
    scalar  k         = x_a[m-2];  // keff
    scalar  lambda    = x_a[m-1];  // current eigenvalue
    scalar  fpfk      = 0.0;
    scalar  fpflambda = 0.0;
    SermentVector Jinc( m-2 );
    SermentVector fJ( m-2 );
    SermentVector fpfJ( m-2 );
    SermentVector temp( m-2 );
    SermentVector MR1J( m-2 );
    SermentVector MR2J( m-2 );
    SermentVector tempC( m );

    VecPlaceArray( Jinc.V, x_a   );
    VecPlaceArray(   fJ.V, f_a   );
    VecPlaceArray( fpfJ.V, fpf_a );

    // we first compute the quantities that need unperturbed responses

  //  cout << " ---begin--- " << endl;
    scalar t = MPI_Wtime();
    problem->R->updateData( k );
    problem->L.updateData( k );
    problem->F.updateData( k );
    problem->A.updateData( k );
  //  cout << " update time for k = " << MPI_Wtime()-t << " seconds " << endl;
   // cout << " ---end--- " << endl;

    // (1)
    // M*R*fJ-lambda*fJ - flambda*Jinc
    problem->R->matVec( fJ,   temp );    // = R*f(1:m-2)
    problem->M->matVec( temp, fpfJ );    // = M*R*f(1:m-2)
    fpfJ.vecAYPV( none*lambda, fJ );    // = M*R*f(1:m-2) - lambda*f(1:m-2)   
    fpfJ.vecAYPV( -flambda, Jinc  );    // = (M*R-lambda*I)*fJ-flambda*J

    // M*R_k*Jinc, pre-perturb
    problem->R->matVec( Jinc, temp );
    problem->M->matVec( temp, MR1J );   // M*R(1)*J

    // (F_k-k*L_k)*Jinc*f(m-1) -- computation by parts:
    scalar F1 =  problem->F.vecDot(Jinc);
    scalar L1 =  problem->A.vecDot(Jinc) + problem->L.computeLeakage(Jinc);

    // (F-k*L)*fJ
    fpfk = problem->F.vecDot( fJ );                         
    fpfk = fpfk -  k*( problem->A.vecDot(fJ) +     // = (F-k*L)*fJ
           problem->L.computeLeakage(fJ) ) ;

    // -Jinc'*fJ
    fpflambda = -Jinc.vecDot(fJ);

//    ///
//    Mat MR;
//    Mat Rcsr;
//    MatConvert(problem->R->M, MATAIJ, MAT_INITIAL_MATRIX, &Rcsr );
//    MatMatMult( problem->M->M, Rcsr, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &MR );
//    SermentVector D(m-2);
//    for (integer i = 0; i < m-2; i++)
//        D.insertVal( i, -lambda );
//    MatDiagonalSet( MR, D.V, ADD_VALUES );
//    //MatView( MR, PETSC_VIEWER_STDOUT_SELF );

//    MatMult( MR, fJ.V, D.V );
    
    // now the perturbed quantities -- update R,L,F,A with new k
    scalar epsilon = 1e-8;  // ~ sqrt(eps_mach), decent static value

   // cout << " ---2begin--- " << endl;
    t = MPI_Wtime();
    problem->R->updateData( k+epsilon );
    problem->L.updateData( k+epsilon );
    problem->F.updateData( k+epsilon );
    problem->A.updateData( k+epsilon );
   // cout << " update time for k = " << MPI_Wtime()-t << " seconds " << endl;
   // cout << " ---2end--- " << endl;

    // jac(1,2), post-perturb
    problem->R->matVec( Jinc, temp );
    problem->M->matVec( temp, MR2J ); // = M*R(2)*Jinc
    MR2J.vecAYPV( none, MR1J );      // = M*(R(2)-R(1))*Jinc
    MR2J.vecScale( 1.0/epsilon );    // ~ M*R_k*Jinc

    fpfJ.vecAYPV( fk, MR2J ); //= M*R-lambda*I*fJ - flambda*Jinc + M*R_k*Jinc*fk

//    // fpfJ =  D.V + -flambda*J + fk*tempB
//    for (integer i = 0; i < m-2; i++)
//        temp.insertVal( i, 0 );
//    temp.vecAYPV( one, D );
//    temp.vecAYPV( fk, tempB );
//    temp.vecAYPV( -flambda, Jinc );

//    cout << " Jinc = " << endl;
//    Jinc.viewMe();
//    cout << " fpfJ in temp = " << endl;
//    temp.viewMe();

    scalar F2 =  problem->F.vecDot(Jinc);
    scalar L2 =  problem->A.vecDot(Jinc) + problem->L.computeLeakage(Jinc);

    fpfk = fpfk + ( (F2-F1)/epsilon - k*(L2-L1)/epsilon - L1 )*fk;

    fpf_a[m-2] = fpfk;
    fpf_a[m-1] = fpflambda;

    // we need to restore x, f, and fpf
    VecResetArray( Jinc.V );
    VecResetArray(   fJ.V );
    VecResetArray( fpfJ.V );
    VecRestoreArray( x->V, &x_a );
    VecRestoreArray( f,    &f_a   );
    VecRestoreArray( fpf,  &fpf_a );

    if (debug)
    {
        cout << " fpf via my jacobian... " << endl;
        VecView( fpf, PETSC_VIEWER_STDOUT_SELF );
        MatMult( JJ, f, tempC.V );
        cout << " fpf via fd jacobian... " << endl;
        tempC.viewMe();
    }

//    cout << " ********** END JACOBIAN ********** " << endl;
//    cout << " x is !! " << endl;    
//    VecView( x->V, PETSC_VIEWER_STDOUT_SELF );
//    cout << " f is !! " << endl;    
//    VecView( f, PETSC_VIEWER_STDOUT_SELF );
//    cout << " fpf is !! " << endl;    
//    VecView( fpf, PETSC_VIEWER_STDOUT_SELF );
//    cout << " ********** END JACOBIAN ********** " << endl;
    Jinc.releaseMe();
    fJ.releaseMe();
    fpfJ.releaseMe();

    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief This is a matrix-vector multiplication wrapper function for use with
 *        JacobianShell.
 *
 * It exists outside the SermentMatrixShell class so that it can be passed
 * to Petsc, as Petsc needs pointers to functions but doesn't like 
 * pointers to member functions.  Here, the context associated with the instance
 * of JacobianShell is retrieved.  That context is simply a pointer to the 
 * global (Newton) solver in disguise, so we cast it into the correct form and
 * call JacobianShell's own member function for multiplication.
 *
 */
PetscErrorCode myMatVecWrap( Mat M, Vec f, Vec fpf )
{

    void *ctx;
    MatShellGetContext( M, &ctx );   // get the context set at initialization
    Newton *me = (Newton*)ctx;       // and cast it as a pointer to Newton

    me->fp->myMatVec( M, f, fpf );   // perform the action

    return 0;
}

//---------------------------------------------------------------------------//
//                 end of JacobianShell.cc
//---------------------------------------------------------------------------//

