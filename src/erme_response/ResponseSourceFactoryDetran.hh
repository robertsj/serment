//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ResponseSourceFactoryDetran.hh
 *  @author robertsj
 *  @date   Aug 31, 2012
 *  @brief  ResponseSourceFactoryDetran build specialization
 */
//---------------------------------------------------------------------------//

#ifndef erme_response_RESPONSESOURCEFACTORYDETRAN_HH_
#define erme_response_RESPONSESOURCEFACTORYDETRAN_HH_

#include "erme_geometry/CartesianNodeDetran.hh"
#include "transport/DimensionTraits.hh"
#include "boundary/BoundaryDiffusion.hh"
#include "boundary/BoundarySN.hh"
#include "boundary/BoundaryMOC.hh"

namespace erme_response
{

//---------------------------------------------------------------------------//
inline ResponseSourceFactory::SP_source
ResponseSourceFactory::build_detran(SP_node node, SP_indexer indexer)
{
  erme_geometry::CartesianNodeDetran::SP_node detran_node = node;
  Require(detran_node->db());
  Require(detran_node->material());
  Require(detran_node->mesh());

  using namespace detran;

  size_t d = detran_node->mesh()->dimension();
  std::string eq = "dd";
  size_t discretization = 0;
  if (detran_node->db()->check("equation"))
    eq = detran_node->db()->get<std::string>("equation");
  if (eq == "scmoc" || eq == "ddmoc")
    discretization = 1;
  else if (eq == "diffusion")
    discretization = 2;

  SP_source s;
  if (d == 1)
  {
    if (discretization == 0)
    {
      s = new ResponseSourceDetran<BoundarySN<_1D> >(node, indexer);
    }
    else if (discretization == 2)
    {
      s = new ResponseSourceDetran<BoundaryDiffusion<_1D> >(node, indexer);
    }
  }
  else if (d == 2)
  {
    if (discretization == 0)
    {
     // s = new ResponseSourceDetran<BoundarySN<_2D> >(node, indexer);
    }
    else if (discretization == 1)
    {
     // s = new ResponseSourceDetran<BoundaryMOC<_2D> >(node, indexer);
    }
    else if (discretization == 2)
    {
      s = new ResponseSourceDetran<BoundaryDiffusion<_2D> >(node, indexer);
    }
  }
  else if (d == 3)
  {
    if (discretization == 0)
    {
     // s = new ResponseSourceDetran<BoundarySN<_3D> >(node, indexer);
    }
    else if (discretization == 2)
    {
      s = new ResponseSourceDetran<BoundaryDiffusion<_3D> >(node, indexer);
    }
  }

  Ensure(s);
  return s;
}

} // end namespace erme_response

#endif /* erme_response_RESPONSESOURCEFACTORYDETRAN_HH_ */
