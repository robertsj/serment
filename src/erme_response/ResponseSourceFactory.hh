//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ResponseSourceFactory.hh
 *  @brief  ResponseSourceFactory
 *  @author Jeremy Roberts
 *  @date   Aug 29, 2012
 */
//---------------------------------------------------------------------------//

#ifndef erme_response_RESPONSESOURCEFACTORY_HH_
#define erme_response_RESPONSESOURCEFACTORY_HH_

#include "ResponseIndexer.hh"
#include "ResponseSource.hh"
// Concrete source types
#include "ResponseSourceDummy.hh"
#include "ResponseSourceDetran.hh"
#include "ResponseSourceDatabase.hh"

namespace erme_response
{

/**
 *  @class ResponseSourceFactory
 *  @brief Constructs response
 *
 */
class ResponseSourceFactory
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef ResponseSource::SP_source     SP_source;
  typedef ResponseSource::SP_node       SP_node;
  typedef ResponseIndexer::SP_indexer   SP_indexer;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /**
   *  @brief Build a response source
   *
   *  This is a factory method that calls an implemention for
   *  each type of node to be constructed.
   *
   *  @param node   Smart pointer to node
   */
  template <typename SP_NODE>
  SP_source build(SP_NODE node, SP_indexer indexer)
  {
    // Preconditions
    Require(node);

    SP_source s;

    // Try to dynamic-cast to all the node types available
    if (dynamic_cast<erme_geometry::CartesianNodeDummy*>(node.bp()))
      s = build_dummy(node, indexer);
    else if(dynamic_cast<erme_geometry::CartesianNodeDetran*>(node.bp()))
      s = build_detran(node, indexer);
    else if(dynamic_cast<erme_geometry::CartesianNode*>(node.bp()))
      s = build_database(node, indexer);
    else
      THROW("Unsupported response source type");

    // Postconditions
    Ensure(s);
    return s;
  }

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  SP_source build_detran(SP_node node, SP_indexer indexer);
  SP_source build_dummy(SP_node node, SP_indexer indexer);
  SP_source build_database(SP_node node, SP_indexer indexer);
  //SP_source build_openmc(SP_node node, SP_indexer indexer);

};

} // end namespace erme_response


//---------------------------------------------------------------------------//
// FACTORY METHOD SPECIALIZATIONS
//---------------------------------------------------------------------------//

#include "ResponseSourceFactoryDummy.hh"
#include "ResponseSourceFactoryDetran.hh"
#include "ResponseSourceFactoryDatabase.hh"

#endif // erme_response_RESPONSESOURCEFACTORY_HH_

//---------------------------------------------------------------------------//
//              end of file ResponseSourceFactory.hh
//---------------------------------------------------------------------------//
