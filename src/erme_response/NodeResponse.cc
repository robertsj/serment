//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  NodeResponse.cc
 *  @brief NodeResponse member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "NodeResponse.hh"
#include <cstdio>

namespace erme_response
{

//----------------------------------------------------------------------------//
NodeResponse::NodeResponse(const size_t N,
                           const size_t number_surfaces,
                           const size_t number_pins)
  : d_N(N)
  , d_number_surfaces(number_surfaces)
  , d_number_pins(number_pins)
  , d_boundary_response(N, vec_dbl(N, 0.0))
  , d_fission_response(N, 0.0)
  , d_absorption_response(N, 0.0)
  , d_leakage_response(N, vec_dbl(number_surfaces, 0.0))
  , d_nodal_power(N, 0.0)
  , d_pin_power(N, vec_dbl(number_pins, 0.0))
{
  Require(d_N > 1); // At least 2 surfaces, each with at least 1 moment
  Require(d_number_surfaces > 1);
}

//----------------------------------------------------------------------------//
void NodeResponse::display() const
{
  printf("\n");
  printf("Boundary Responses \n");
  printf("------------------ \n");
  for (int out = 0; out < d_N; out++)
  {
    for (int in = 0; in < d_N; in++)
    {
      printf("%16.8e ", d_boundary_response[in][out]);
    }
    printf("\n");
  }
  printf("\n");
  printf("Leakage Responses \n");
  printf("------------------ \n");
  for (int out = 0; out < d_number_surfaces; out++)
  {
    for (int in = 0; in < d_N; in++)
    {
      printf("%16.8e ", d_leakage_response[in][out]);
    }
    printf("\n");
  }
  printf("\n");
  printf("Fission Responses \n");
  printf("------------------ \n");
  for (int in = 0; in < d_N; in++)
  {
    printf("%16.8e ", d_fission_response[in]);
  }
  printf("\n\n");
  printf("Absorption Responses \n");
  printf("------------------ \n");
  for (int in = 0; in < d_N; in++)
  {
    printf("%16.8e ", d_absorption_response[in]);
  }
  printf("\n\n");
  printf("Nodal Power Responses \n");
  printf("------------------ \n");
  for (int in = 0; in < d_N; in++)
  {
    printf("%16.8e ", d_nodal_power[in]);
  }
  printf("\n\n");
  printf("Pin Power Responses \n");
  for (int p = 0; p < d_number_pins; p++)
  {
    for (int in = 0; in < d_N; in++)
    {
      printf("%16.8e ", d_pin_power[in][p]);
    }
    printf("\n");
  }
  printf("\n");
}

} // end namespace erme_response

//----------------------------------------------------------------------------//
//              end of file NodeResponse.cc
//----------------------------------------------------------------------------//
