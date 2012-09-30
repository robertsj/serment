//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Connect.hh
 * \brief  Connect class definition
 * \author Jeremy Roberts
 * \date   Aug 23, 2012
 */
//---------------------------------------------------------------------------//

#ifndef CONNECT_HH_
#define CONNECT_HH_

#include "linear_algebra/Matrix.hh"
#include "erme_geometry/NodeList.hh"
#include "erme_response/ResponseIndexer.hh"

namespace erme
{

/*!
 *  \class Connect
 *  \brief Defines the geometric relationship between nodes
 *
 *  In the simplest case of one unknown per node surface, the connectivity
 *  matrix is simply an adjacency matrix in which a 1 is placed at any
 *  row and column location representing the surface of one node adjoining
 *  the surface of another node.
 *
 *  Vacuum or reflective boundary conditions are built into the connectivity
 *  matrix.  Reflective conditions require a polarity shift for odd
 *  moments; this applies to some combinations of space and angle basis
 *  functions.
 */
class Connect: public linear_algebra::Matrix
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran::SP<Connect>                         SP_connect;
  typedef erme_geometry::NodeList::SP_nodelist        SP_nodelist;
  typedef erme_response::ResponseIndexer::SP_indexer  SP_indexer;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   *  \param nodes    Pointer to list of nodes
   *  \param indexer  Pointer to response indexer
   */
  Connect(SP_nodelist nodes, SP_indexer indexer);

private:

  //-------------------------------------------------------------------------//
  // PRIVATE DATA
  //-------------------------------------------------------------------------//

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

};


} // end namespace detran

#endif // CONNECT_HH_ 

//---------------------------------------------------------------------------//
//              end of file Connect.hh
//---------------------------------------------------------------------------//
