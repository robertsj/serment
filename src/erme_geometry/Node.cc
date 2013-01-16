//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Node.cc
 * \brief  Node 
 * \author Jeremy Roberts
 * \date   Aug 22, 2012
 */
//---------------------------------------------------------------------------//

#include "Node.hh"
#include "NodeSerialization.hh"

namespace erme_geometry
{

//---------------------------------------------------------------------------//
Node::Node(const size_t  dimension,
           const size_t  number_surfaces,
           std::string   name,
           vec2_size_t   so,
           vec_size_t    po,
           vec_size_t    ao,
           vec_size_t    eo)
  : d_dimension(dimension)
  , d_number_surfaces(number_surfaces)
  , d_name(name)
  , d_spatial_order(so)
  , d_polar_order(po)
  , d_azimuthal_order(ao)
  , d_energy_order(eo)
{
  // Preconditions
  Require(d_dimension > 0 and d_dimension <= 3);
  Require(d_number_surfaces > 0);
  Require(d_spatial_order.size()    == d_number_surfaces);
  Require(d_polar_order.size()      == d_number_surfaces);
  Require(d_azimuthal_order.size()  == d_number_surfaces);
  Require(d_energy_order.size()     == d_number_surfaces);
}

//---------------------------------------------------------------------------//
Node::~Node()
{
  /* ... */
}


void Node::display() const
{
  std::cout << "----------------------------" << std::endl;
  std::cout << "NODE" << std::endl;
  std::cout << "----------------------------" << std::endl;
  std::cout << "  dimension = " << d_dimension << std::endl;
  std::cout << "   surfaces = " << d_number_surfaces << std::endl;
  std::cout << "      name  = " << d_name << std::endl;
  std::cout << "----------------------------" << std::endl;
  for (size_t s = 0; s < d_number_surfaces; ++s)
  {
    std::cout << "  SURFACE " << s << std::endl;
    if (d_dimension > 1)
      std::cout << "    spatial 0 = " << d_spatial_order[s][0] << std::endl;
    if (d_dimension > 2)
      std::cout << "    spatial 1 = " << d_spatial_order[s][1] << std::endl;
    std::cout << "        polar = " << d_polar_order[s] << std::endl;
    std::cout << "    azimuthal = " << d_azimuthal_order[s] << std::endl;
  }
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file Node.cc
//---------------------------------------------------------------------------//
