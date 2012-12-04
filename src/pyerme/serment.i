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
  
// Response
#include "NodeResponse.hh"
#include "ResponseIndex.hh"
#include "ResponseIndexer.hh"

// Linear algebra
#include "LinearAlgebraSetup.hh"
  
// Manager
#include "erme_utils/ManagerERME.hh"

%} // end module pyserment

%include "serment_config.h"

// STL
%include "std_string.i"
%include "std_vector.i"

// Detran
%include "utilities/detran_utilities.i"
%include "material/detran_material.i"
%include "geometry/detran_geometry.i"

// Geometry
%include "erme_geometry/NeighborSurface.hh"
%template(vec_neighbor) std::vector<erme_geometry::NeighborSurface>;
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

//// Response
//%include "NodeResponse.hh"
//%template(NodeResponseSP)  detran_utilities::SP<erme_response::NodeResponse>;
//%include "ResponseIndex.hh"
//%template(ResponseIndexSP)  detran_utilities::SP<erme_response::ResponseIndex>;
//%include "ResponseIndexer.hh"
//%template(ResponseIndexerSP)  detran_utilities::SP<erme_response::ResponseIndexer>;
//
//// PETSc/SLEPc initialization
//%include "linear_algebra/LinearAlgebraSetup.hh"
//
//// Manager
//%include "erme_utils/ManagerERME.hh"
//%template(ManagerERMESP)          detran_utilities::SP<erme_utils::ManagerERME>;
