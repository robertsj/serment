//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   NodeSerialization.hh
 * \brief  Convenience header for holding Boost serialization macros
 * \author Jeremy Roberts
 * \date   Oct 4, 2012
 */
//---------------------------------------------------------------------------//

#ifndef NODESERIALIZATION_HH_
#define NODESERIALIZATION_HH_

#include "CartesianNodeDetran.hh"
#include "DummyNode.hh"

// This header must be included by node list

// Node
BOOST_SERIALIZATION_ASSUME_ABSTRACT(erme_geometry::Node)
// Cartesian
BOOST_CLASS_EXPORT_KEY(erme_geometry::CartesianNode)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(erme_geometry::CartesianNode)
// Hexagonal
//BOOST_CLASS_EXPORT_KEY(erme_geometry::HexagonalNode)
//BOOST_SERIALIZATION_ASSUME_ABSTRACT(erme_geometry::HexagonalNode)
// Detran specialization
BOOST_CLASS_EXPORT_KEY(erme_geometry::CartesianNodeDetran)
// Dummy specialization
BOOST_CLASS_EXPORT_KEY(erme_geometry::CartesianNodeDummy)


#endif // NODESERIALIZATION_HH_ 

//---------------------------------------------------------------------------//
//              end of file NodeSerialization.hh
//---------------------------------------------------------------------------//
