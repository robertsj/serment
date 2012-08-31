//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ResponseSourceFactoryDetran.hh
 * \author robertsj
 * \date   Aug 31, 2012
 * \brief  ResponseSourceFactoryDetran build specialization
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef RESPONSESOURCEFACTORYDETRAN_HH_
#define RESPONSESOURCEFACTORYDETRAN_HH_

#include "erme_geometry/CartesianNodeDetran.hh"
#include "ResponseSourceDetran.hh"

namespace erme_response
{

template <>
ResponseSourceFactory::SP_source
ResponseSourceFactory::build(detran::SP<erme_geometry::CartesianNodeDetran> node)
{
  SP_source s;
  return s;
}

} // end namespace erme_response


#endif /* RESPONSESOURCEFACTORYDETRAN_HH_ */
