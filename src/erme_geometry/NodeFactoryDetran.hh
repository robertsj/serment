//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  NodeFactoryDetran.hh
 *  @brief NodeFactoryDetran class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef erme_geometry_NODEFACTORYDETRAN_HH_
#define erme_geometry_NODEFACTORYDETRAN_HH_

#include "NodeFactory.hh"

namespace erme_geometry
{

/**
 *  @class NodeFactoryDetran
 *  @brief Build Detran-based nodes
 */
class NodeFactoryDetran: public NodeFactory
{

public:

  /// Constructor
  NodeFactoryDetran(){}

  /// Virtual destructor
  virtual ~NodeFactoryDetran(){}

  /**
   *  @brief Create a node
   *  @param db   Parameter database.
   */
  SP_node create_node(SP_db db, SP_material material, SP_mesh mesh);

private:

};

} // end namespace erme_geometry

#endif // erme_geometry_NODEFACTORYDETRAN_HH_

//----------------------------------------------------------------------------//
//              end of file NodeFactoryDetran.hh
//----------------------------------------------------------------------------//
