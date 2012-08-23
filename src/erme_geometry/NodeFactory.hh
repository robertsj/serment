//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   NodeFactory.hh
 * \brief  NodeFactory class definition
 * \author Jeremy Roberts
 * \date   Aug 23, 2012
 */
//---------------------------------------------------------------------------//

#ifndef NODEFACTORY_HH_
#define NODEFACTORY_HH_

#include "InputDB.hh"
#include "Node.hh"

namespace erme_geometry
{

/*!
 *  \class NodeFactory
 *  \brief Base class for constructing various node types
 */
class NodeFactory
{

public:

  typedef detran::InputDB::SP_input SP_db;
  typedef Node::SP_node SP_node;

  /// Constructor
  NodeFactory(){}

  /// Virtual destructor
  virtual ~NodeFactory(){}

  /*!
   *  \brief Create a node
   *
   *  Everything needed to build the node must
   *  exit in the database.
   *
   *  \param db   Parameter database.
   */
  virtual SP_node create_node(SP_db db) = 0;

private:

};

} // end namespace erme_geometry

#endif // NODEFACTORY_HH_ 

//---------------------------------------------------------------------------//
//              end of file NodeFactory.hh
//---------------------------------------------------------------------------//
