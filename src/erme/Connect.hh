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
 *  Boundary conditions are built into
 */
class Connect: public linear_algebra::Matrix
{

public:

  typedef erme_geometry::NodeList         NodeList;
  typedef erme_response::ResponseIndexer  ResponseIndexer;

  Connect(erme_geometry::NodeList &nodes,
          erme_response::ResponseIndexer &indexer);

private:
//
//  const NodeList&         d_nodes;
//  const ResponseIndexer&  d_indexer;

};


} // end namespace detran

#endif // CONNECT_HH_ 

//---------------------------------------------------------------------------//
//              end of file Connect.hh
//---------------------------------------------------------------------------//
