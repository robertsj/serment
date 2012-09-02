//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseSourceFactoryDummy.hh
 * \author robertsj
 * \date   Aug 31, 2012
 * \brief  ResponseSourceFactory build specialization for Dummy nodes
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef RESPONSESOURCEFACTORYDUMMY_HH_
#define RESPONSESOURCEFACTORYDUMMY_HH_

#include "erme_geometry/DummyNode.hh"
#include "ResponseSourceDummy.hh"

namespace erme_response
{

inline ResponseSourceFactory::SP_source
ResponseSourceFactory::build_dummy(SP_node node)
{
  SP_source s(new ResponseSourceDummy(node));
  return s;
}

} // end namespace erme_response

#endif /* RESPONSESOURCEFACTORYDUMMY_HH_ */
