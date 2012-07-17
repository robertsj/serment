//----------------------------------*-C++-*----------------------------------//
/*
 * \file   GlobalInput.hh
 * \author Jeremy Roberts
 * \date   10/26/2010
 * \brief  GlobalInput class definition.
 * \note   Copyright (C) 2010 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 152                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-09-27 06:34:56 -0400 (Tue, 27 Sep 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef GLOBALINPUT_HH
#define GLOBALINPUT_HH

#include <iostream>
#include <fstream>
#include <vector>
#include "petscvec.h"
#include "petscsys.h"
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>
#include <libxml/relaxng.h>
#include "LinAlg.hh"
#include "utilities/InputXML.hh"
#include "utilities/SP.hh"

using namespace std;

//===========================================================================//
/*
 * \class GlobalInput
 * \brief This class handles all input for Serment.
 *
 *  This class used the libxml2 library to parse xml inputs.  Some of the code
 *  is very much inspired by K. Huff's cyclus input processing routines.
 *
 * \todo Finish the map-based input utility.
 *
 */
//===========================================================================//
class GlobalInput : public InputXML
{
public:

  /// Typedefs
  //\{
  typedef typename util::SP<GlobalInput>    SP_globalinput;
  //\}

  /// ------------------------------------------------------------------
  /// GLOBAL CONTROL
  /// ------------------------------------------------------------------

  /// \brief problem input file name
  string file;

  /// \brief problem description
  string disc;

  /*! \brief response function source
   *
   * <UL>
   *    <LI> rfsource = 0 = database </LI>
   *    <LI> rfsource = 1 = diff2d </LI>
   * </UL>
   */
  integer rfsource;

  /*! \brief response function sink
   *
   * <UL>
   *    <LI> rfsink = 0 = database </LI>
   *    <LI> rfsink = 1 = power iteration </LI>
   *    <LI> rfsink = 2 = power iteration + aitken </LI>
   *    <LI> rfsink = 3 = newton-krylov </LI>
   *    <LI> rfsink = 4 = jacobian-free newton-krylov </LI>
   * </UL>
   *
   */
  integer rfsink;

  /// \brief source file name (or dbname for generating)
  string rfsourcefile;

  /// \brief database file for output
  string dbfile;

  /// \brief \f$ k_{e\!f\!f} \f$ iteration relative precision (only ptype = 1)
  scalar epsk;

  /// \brief flux iteration relative precision (or source precision, ptype = 1)
  scalar epss;

  /// \brief maximum number of outer iterations
  scalar maxit;

  /// \brief initial guess for \f$ k_{e\!f\!f} \f$
  scalar keff;

  /// \brief number of element types
  integer numtypes;

  /// \brief number of elements
  integer numel;

  /// \brief energy groups
  integer numgroups;

  /// \brief spatial expansion order
  integer spaceord;

  /// \brief angular expansion order
  integer angleord;

  /// \brief number of faces per element
  integer faces;

  /// \brief degrees of freedom
  integer degfree;

  /// \brief preconditioner type (0=none,1=ilu)
  integer pctype;

  /// \brief ilu factor level
  integer ilulevel;

  // ------------------------------------------------------------------
  // OUTPUT (NOT FULLY KNOWN YET)
  // ------------------------------------------------------------------

  /// \brief should I print? 1=yes, 0=no
  integer printout;

  /// \brief where should I print?
  string outfile;

  /// \brief should i plot? 1=yes, 0=no
  integer plotout;

  /// \brief where should I print?
  string fluxfile;

  // ------------------------------------------------------------------
  // GLOBAL GEOMETRY
  // ------------------------------------------------------------------

  /// \brief left (x) global boundary condition
  integer bcl;

  /// \brief right (x)  global boundary condition
  integer bcr;

  /// \brief bottom (y) global boundary condition
  integer bcb;

  /// \brief top (y) global boundary condition
  integer bct;

  /// \brief north (z)
  integer bcn;

  /// \brief south (z)
  integer bcs;

  /// \brief number of elements in x direction
  integer elemx;

  /// \brief number of elements in y direction
  integer elemy;

  /// \brief number of elements in z direction
  integer elemz;

  /// \brief element placement [][][]
  vector<vector<vector<scalar> > > elements;

  GlobalInput(string S);

  ~GlobalInput();

  void echoInput();

private:
  // should be more private things...

  bool processFile();

  void sayHello();

};

#endif // GLOBALINPUT_HH
//---------------------------------------------------------------------------//
//                 end of GlobalInput.hh
//---------------------------------------------------------------------------//

