//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseSourceFactory.hh
 * \brief  ResponseSourceFactory 
 * \author Jeremy Roberts
 * \date   Aug 29, 2012
 */
//---------------------------------------------------------------------------//

#ifndef RESPONSESOURCEFACTORY_HH_
#define RESPONSESOURCEFACTORY_HH_

#include "ResponseSource.hh"

namespace erme_response
{

/*!
 *  \class ResponseSourceFactory
 *  \brief Constructs response
 *
 */
class ResponseSourceFactory
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef ResponseSource::SP_source SP_source;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /// Build a node.
  template <typename NODE>
  SP_source build(detran::SP<NODE> node)
  {
    THROW("The source factory build method must be specialized.");
  }

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

};

} // end namespace erme_response


//---------------------------------------------------------------------------//
// FACTORY METHOD SPECIALIZATIONS
//---------------------------------------------------------------------------//

#include "ResponseSourceFactoryDummy.hh"
//#include "ResponseSourceFactoryDetran.hh"
//#include "ResponseSourceFactoryOpenMC.hh"

#endif // RESPONSESOURCEFACTORY_HH_ 

//---------------------------------------------------------------------------//
//              end of file ResponseSourceFactory.hh
//---------------------------------------------------------------------------//
