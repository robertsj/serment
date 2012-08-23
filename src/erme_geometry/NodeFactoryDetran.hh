//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   NodeFactoryDetran.hh
 * \brief  NodeFactoryDetran 
 * \author Jeremy Roberts
 * \date   Aug 23, 2012
 */
//---------------------------------------------------------------------------//

#ifndef NODEFACTORYDETRAN_HH_
#define NODEFACTORYDETRAN_HH_

#include "NodeFactory.hh"

namespace erme_geometry
{

/*!
 *  \class NodeFactoryDetran
 *  \brief Build Detran-based nodes
 */
class NodeFactoryDetran: public NodeFactory
{

public:

  /// Constructor
  NodeFactoryDetran(){}

  /// Virtual destructor
  virtual ~NodeFactoryDetran(){}

  /*!
   *  \brief Create a node
   *  \param db   Parameter database.
   */
  SP_node create_node(SP_db db,
                      SP_material material,
                      SP_mesh mesh);

private:

};

} // end namespace erme_geometry

#endif // NODEFACTORYDETRAN_HH_ 

//---------------------------------------------------------------------------//
//              end of file NodeFactoryDetran.hh
//---------------------------------------------------------------------------//
