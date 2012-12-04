//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Node.i.hh
 *  @author robertsj
 *  @date   Oct 1, 2012
 *  @brief  Node.i class definition.
 */
//---------------------------------------------------------------------------//

#ifndef NODE_I_HH_
#define NODE_I_HH_

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
  // Preconditions
  Require(s < d_number_surfaces);
  Require(d < d_spatial_order[s].size());
  return d_spatial_order[s][d];
}

//---------------------------------------------------------------------------//
inline Node::size_t Node::polar_order(const size_t s) const
{
  // Preconditions
  Require(s < d_number_surfaces);
  return d_polar_order[s];
}

//---------------------------------------------------------------------------//
inline Node::size_t Node::azimuthal_order(const size_t s) const
{
  // Preconditions
  Require(s < d_number_surfaces);
  return d_azimuthal_order[s];
}

//---------------------------------------------------------------------------//
inline Node::size_t Node::energy_order(const size_t s) const
{
  // Preconditions
  Require(s < d_number_surfaces);
  return d_energy_order[s];
}

} // end namespace erme_geometry

#endif /* NODE_I_HH_ */
