//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Node.i.hh
 *  @author robertsj
 *  @date   Oct 1, 2012
 *  @brief  Node.i class definition.
 */
//---------------------------------------------------------------------------//

#ifndef erme_geometry_NODE_I_HH_
#define erme_geometry_NODE_I_HH_

namespace erme_geometry
{

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
inline Node::size_t Node::dimension() const
{
  return d_dimension;
}

//---------------------------------------------------------------------------//
inline Node::size_t Node::number_surfaces() const
{
  return d_number_surfaces;
}

//---------------------------------------------------------------------------//
inline std::string Node::name() const
{
  return d_name;
}

//---------------------------------------------------------------------------//
inline Node::size_t Node::spatial_order(const size_t s, const size_t d) const
{
  Require(s < d_number_surfaces);
  size_t v = 0;
  if (d_spatial_order.size())
  {
    Require(d < d_spatial_order[s].size());
    v = d_spatial_order[s][d];
  }
  return v;
}

//---------------------------------------------------------------------------//
inline Node::size_t Node::polar_order(const size_t s) const
{
  Require(s < d_number_surfaces);
  return d_polar_order[s];
}

//---------------------------------------------------------------------------//
inline Node::size_t Node::azimuthal_order(const size_t s) const
{
  Require(s < d_number_surfaces);
  return d_azimuthal_order[s];
}

//---------------------------------------------------------------------------//
inline Node::size_t Node::energy_order(const size_t s) const
{
  Require(s < d_number_surfaces);
  return d_energy_order[s];
}

} // end namespace erme_geometry

#endif /* erme_geometry_NODE_I_HH_ */
