//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   TwoGroupKernelSolver.hh
 * \author Jeremy Roberts
 * \date   10/14/2011
 * \brief  TwoGroupKernelSolver class definition.
 * \note   Copyright (C) 2011 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 140                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-09-14 12:53:40 -0400 (Wed, 14 Sep 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#include "petscvec.h"
#include "petscmat.h"
#include "petscsys.h"
#include "petscksp.h"
#include "petscpc.h"

//#include "ResponseFunction.hh"
#include "utilities/LegendrePoly.hh"
#include "utilities/SP.hh"

namespace erme
{

//===========================================================================//
/*!
 * \class TwoGroupKernel
 * \brief Simple diff2d alternative that borrows some of its code.
 *
 * This aims to be a simple, fast solver limited to two group, homogeneous
 * nodes without upscatter.
 *
 */
//===========================================================================//

class TwoGroupKernelSolver
{

public:

  /// \name Typedefs
  //\{
  typedef util::SP<TwoGroupKernelSolver>    SP_kernelsolver;
  //typedef util::SP<ResponseFunction>        SP_responsefunction;
  typedef std::vector<double>               Vec_Dbl;
  typedef std::vector<int>                  Vec_Int;
  typedef std::vector<Vec_Dbl>              Mat_Dbl;
  //\}

  /*!
   *  \brief Constructor
   *
   *  Initializes the matrix, solvers, etc.
   *
   *  \param    nh      number of spatial mesh
   *  \param    order   order of responses to generate
   *
   */
  TwoGroupKernelSolver(int nh, int order);

  /*!
   *  \brief Update the matrix
   *
   *  \param    data    group data for this state
   */
  void update_matrix(Vec_Dbl &data);


  void set_source(Vec_Dbl &data, int group, int order);
  /*!
   *  \brief Solve the system
   *
   */
  void solve();


private:

  /// \name Local Response Storage
  //\{

  /// Current response.
  Vec_Dbl d_current_rf;

  /// Leakage response.
  Vec_Dbl d_leakage_rf;

  /// Fission response.
  Vec_Dbl d_fission_rf;

  /// Absorption response.
  Vec_Dbl d_absorption_rf;

  /// One group flux response.
  Vec_Dbl d_total_flux_rf;

  /// Group flux response.
  Vec_Dbl d_group_flux_rf;

  /// Group flux gradient response.
  Vec_Dbl d_gradient_rf;

  //\}

  /// \name System Data
  //\{

  /// Number spatial divisions.
  int d_nh;

  /// Matrix size.
  int d_n;

  /// Number of nonzeros in the matrix.
  int d_nnz;

  /// row pointers
  Vec_Int d_rowindex;

  /// row count
  Vec_Int d_rowcount;

  /// column index
  Vec_Int d_column;


  //}

  /// \name PETSc matrices and solvers.
  //\{

  // System matrix.
  Mat A;

  // Preconditioner.
  PC prec;

  // Solver context.
  KSP ksp;

  // Vectors  and Dummy arrays
  Vec v_f;          double *a_f;
  Vec v_si;         double *a_si;
  Vec v_sf;         double *a_sf;
  Vec v_b;          double *a_b;
  Vec v_phi;        double *a_phi;
  Vec v_phi_old;    double *a_phi_old;
  Vec v_err;        double *a_err;
  double *a_A; // matrix dummy

  LegendrePoly d_legendre;

  //}

  /// \name Implementation
  //\{
  void matrix(Vec_Dbl &data, double *t, double *f);

  /*!
   *  Compute the matrix column indices.
   *
   */
  void columns(Vec_Int c, Vec_Int rp);

  /*!
   *  Compute the matrix count per row.
   *
   */
  int countrow(int row);

  /*!
   *  Compute the matrix row indices
   *
   */
  void getrowcount();
  void getrowindex();


  //\}

};




} // end namespace erme

//---------------------------------------------------------------------------//
//                 end of TwoGroupKernel.hh
//---------------------------------------------------------------------------//

