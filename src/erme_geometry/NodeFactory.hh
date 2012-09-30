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

#include "Node.hh"
#include "utilities/InputDB.hh"
#include "material/detran_material.hh"
#include "geometry/detran_geometry.hh"

namespace erme_geometry
{

/*!
 *  \class NodeFactory
 *  \brief Base class for constructing various node types
 */
class NodeFactory
{

public:

  /// Useful typedefs
  typedef detran_utilities::InputDB::SP_input       SP_db;
  typedef detran_material::Material::SP_material    SP_material;
  typedef detran_geometry::Mesh::SP_mesh            SP_mesh;
  typedef Node::SP_node                             SP_node;

  /// Constructor
  NodeFactory(){}

  /// Virtual destructor
  virtual ~NodeFactory(){}

  /*!
   *  \brief Create a node
   *
   *  Everything needed to build the node must
   *  exist in the database.
   *
   *  \param db   Parameter database.
   *  \param db
   */
  virtual SP_node create_node(SP_db db,
                              SP_material material = SP_material(),
                              SP_mesh mesh = SP_mesh()) = 0;

private:

};

} // end namespace erme_geometry

#endif // NODEFACTORY_HH_ 

//---------------------------------------------------------------------------//
//              end of file NodeFactory.hh
//---------------------------------------------------------------------------//
