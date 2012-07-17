//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Diff2dInput.hh
 * \author Jeremy Roberts
 * \date   10/26/2010
 * \brief  A class for handling input for diff2d problems.
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 171                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-11-19 20:36:52 -0500 (Sat, 19 Nov 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef DIFF2DINPUT_HH
#define DIFF2DINPUT_HH

#include <iostream>
#include <fstream>
#include <vector>
#include "petscvec.h"
#include "petscsys.h"
#include "../linalg/typedefs.hh"
#include "Diff2dElement.hh"
#include "utilities/SP.hh"
using namespace std;

//===========================================================================//
/*!
 * \class Diff2dInput
 * \brief A class for handling and containing data for diff2d problems.
 *
 *  Diff2dInput is a class that performs functions necessary for both stand-
 *  alone problems using diff2d, a simple 2-D diffusion code, and problems
 *  where diff2d is used to generate response functions.
 *
 */
//===========================================================================//

class Diff2dInput
{
public:

  typedef util::SP<Diff2dInput> SP_input;
  // ------------------------------------------------------------------
  // GLOBAL CONTROL
  // ------------------------------------------------------------------

  /*! \brief problem input file name
   */
  string file;
  /*! \brief problem description
   */
  string disc;
  /*! \brief problem type
   *
   *  ptype = 0 = fixed source (w/ multiplication) \\
         *  ptype = 1 = eigenvalue \\
         *  ptype = 2 = response  \\
         */
  integer ptype;
  /*! \brief k_{e\!f\!f} iteration relative precision (only ptype = 1)
   */
  scalar epsk;
  /*! \brief flux iteration relative precision
   /          (or source precision, ptype = 1)
   */
  scalar epss;
  /*! \brief maximum number of outer iterations
   */
  scalar maxit;
  /*! \brief number of elements (>1 only for ptype = 2)
   */
  integer numel;
  /*! \brief max Legendre order for ptype = 2
   */
  integer maxOrder;

  // ------------------------------------------------------------------
  // MATERIAL DATA
  // ------------------------------------------------------------------

  /// \brief number of energy groups for given data
  integer numg;

  /// \brief number of materials for given data
  integer numm;

  /// \brief material/group diffusion coefficient [][]
  vector<vector<scalar> > dc;

  /// \brief material/group removal cross-section [][]
  vector<vector<scalar> > sr;

  /// \brief material/group absorption cross-section [][]
  vector<vector<scalar> > ab;

  /// \brief material/group $\nu \Sigma_f$ [][]
  vector<vector<scalar> > ns;

  /// \brief material/group chi-spectrum [][]
  vector<vector<scalar> > xi;

  /// \brief material/group scattering cross-section [][][]
  vector<vector<vector<scalar> > > sc;

  /// ------------------------------------------------------------------
  /// OUTPUT (NOT FULLY KNOWN YET)
  /// ------------------------------------------------------------------

  /*! \brief should I print? 1=yes, 0=no
   */
  integer printout;
  /*! \brief where should I print?
   */
  string outfile;
  /*! \brief should i plot? 1=yes, 0=no
   */
  integer plotout;
  /*! \brief where should I print?
   */
  string fluxfile;

  bool debug;

  // -------------------------------------------------------------------
  /*! \brief elements of the problem
   */
  vector<Diff2dElement> elements;

  void readInput(char *file);
  void skipStuff(ifstream &in);

};

#endif // DIFF2DINPUT_HH
//---------------------------------------------------------------------------//
//                 end of Diff2dInput.hh
//---------------------------------------------------------------------------//

