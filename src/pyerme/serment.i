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
#include "erme_response/NodeResponse.hh"
#include "erme_response/ResponseIndex.hh"
#include "erme_response/ResponseIndexer.hh"
// Linear algebra
#include "linear_algebra/LinearAlgebraSetup.hh"
// Manager
#include "erme_solver/ManagerERME.hh"
#include "erme_utils/PostProcess.hh"
#include "erme_utils/Archive.hh"
%} // end module pyserment

%include "output.i"

%include "serment_config.h"

// STL
%include "std_string.i"
%include "std_vector.i"

// Detran
%import "utilities/detran_utilities.i"
%import "angle/detran_angle.i"
%import "material/detran_material.i"
%import "geometry/detran_geometry.i"

// Geometry
%include "erme_geometry/NeighborSurface.hh"
%include "erme_geometry/Node.hh"
%include "erme_geometry/CartesianNode.hh"
%include "erme_geometry/DummyNode.hh"
%include "erme_geometry/CartesianNodeDetran.hh"
%include "erme_geometry/NodeList.hh"
%include "erme_geometry/NodePartitioner.hh"
// Smart pointer templates
%template(NodeListSP)             detran_utilities::SP<erme_geometry::NodeList>;
%template(NodeSP)                 detran_utilities::SP<erme_geometry::Node>;
%template(CartesianNodeSP)        detran_utilities::SP<erme_geometry::CartesianNode>;
%template(CartesianNodeDummySP)   detran_utilities::SP<erme_geometry::CartesianNodeDummy>;
%template(CartesianNodeDetranSP)  detran_utilities::SP<erme_geometry::CartesianNodeDetran>;
// Vector templates
%template(vec_neighbor)  std::vector<erme_geometry::NeighborSurface>;
%template(vec2_neighbor) std::vector<std::vector<erme_geometry::NeighborSurface> >;
%template(vec_point)     std::vector<detran_utilities::Point>;

%inline
{
// Case base node as CartesianNode
detran_utilities::SP<erme_geometry::CartesianNode> 
as_cartesian_node(detran_utilities::SP<erme_geometry::Node> *n)
{
  return detran_utilities::SP<erme_geometry::CartesianNode>(*n);
}
} // end inline


// Response
%include "NodeResponse.hh"
%template(NodeResponseSP)  detran_utilities::SP<erme_response::NodeResponse>;
%include "ResponseIndex.hh"
%template(ResponseIndexSP)  detran_utilities::SP<erme_response::ResponseIndex>;
%include "ResponseIndexer.hh"
%template(ResponseIndexerSP)  detran_utilities::SP<erme_response::ResponseIndexer>;

// PETSc/SLEPc initialization
%include "linear_algebra/LinearAlgebraSetup.hh"

%include "erme_solver/ManagerERME.hh"
%template(ManagerERMESP)          detran_utilities::SP<erme_solver::ManagerERME>;

// ERME utilities
%include "erme_utils/PostProcess.hh"
%template(PostProcessSP)          detran_utilities::SP<erme_utils::PostProcess>;
%include "erme_utils/Archive.hh"
