//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   serment.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for detran library.
 */
//---------------------------------------------------------------------------//

%module pyserment
%{

// Detran core
#include "utilities/detran_utilities.hh"
#include "material/detran_material.hh"
#include "geometry/detran_geometry.hh"

// Geometry
#include "erme_geometry/NeighborSurface.hh"
#include "erme_geometry/DummyNode.hh"
#include "erme_geometry/CartesianNode.hh"
#include "erme_geometry/CartesianNodeDetran.hh"
#include "erme_geometry/NodeFactory.hh"
#include "erme_geometry/NodeFactoryDetran.hh"
#include "erme_geometry/NodeList.hh"

// Manager
#include "erme_utils/Manager.hh"

%} // end module pyserment

%include "serment_config.h"

// STL
%include "std_string.i"
%include "std_vector.i"

// Detran
%include "detran_utilities.i"
%include "detran_geometry.i"
%include "detran_materials.i"


// Geometry
%include "erme_geometry/NeighborSurface.hh"
%include "erme_geometry/Node.hh"
%include "erme_geometry/CartesianNode.hh"
%include "erme_geometry/DummyNode.hh"
%include "erme_geometry/CartesianNodeDetran.hh"
%include "erme_geometry/NodeFactory.hh"
%include "erme_geometry/NodeFactoryDetran.hh"
%include "erme_geometry/NodeList.hh"
//// Smart pointer templates
%template(NodeSP)                 detran_utilities::SP<erme_geometry::Node>;
%template(CartesianNodeSP)        detran_utilities::SP<erme_geometry::CartesianNode>;
%template(CartesianNodeDummySP)   detran_utilities::SP<erme_geometry::CartesianNodeDummy>;
%template(CartesianNodeDetranSP)  detran_utilities::SP<erme_geometry::CartesianNodeDetran>;
%template(NodeListSP)             detran_utilities::SP<erme_geometry::NodeList>;

// Manager
//%include "erme_utils/Manager.hh"
