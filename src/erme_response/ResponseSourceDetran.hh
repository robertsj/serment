//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ResponseSourceDetran.hh
 *  @brief ResponseSourceDetran class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef erme_response_RESPONSESOURCEDETRAN_HH_
#define erme_response_RESPONSESOURCEDETRAN_HH_

#include "erme_response/ResponseSource.hh"
#include "erme_geometry/CartesianNodeDetran.hh"
#include "orthog/OrthogonalBasis.hh"
#include "solvers/FixedSourceManager.hh"
#include "boundary/BoundaryDiffusion.hh"
#include "boundary/BoundarySN.hh"
#include "boundary/BoundaryMOC.hh"
#include "boundary/BoundaryBase.hh"
#include "boundary/BoundaryTraits.hh"
#include "angle/ProductQuadrature.hh"

namespace erme_response
{

/**
 *  @class ResponseSourceDetran
 *  @brief Compute responses using Detran
 *  @tparam B   Boundary type
 */
template <class B>
class ResponseSourceDetran: public ResponseSource
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef erme_geometry::CartesianNodeDetran::SP_node   SP_node;
  typedef detran_utilities::InputDB::SP_input           SP_db;
  typedef detran_material::Material::SP_material        SP_material;
  typedef detran_geometry::Mesh::SP_mesh                SP_mesh;
  typedef detran_orthog::OrthogonalBasis::SP_basis      SP_basis;
  typedef std::vector<SP_basis>                         vec_basis;
  typedef std::vector<vec_basis>                        vec2_basis;
  typedef B                                             Boundary_T;
  typedef typename Boundary_T::SP_boundary              SP_boundary;
  typedef typename B::D_T                               D;
  typedef detran::FixedSourceManager<D>                 Solver_T;
  typedef detran::BoundaryTraits<D>                     BoundaryTraits_T;
  typedef detran::BoundaryValue<D>                      BoundaryValue_T;
  typedef typename detran_utilities::SP<Solver_T>       SP_solver;
  typedef typename Solver_T::SP_quadrature              SP_quadrature;
  typedef detran_angle::ProductQuadrature               ProductQuadrature;
  typedef detran_utilities::SP<ProductQuadrature>       SP_productquadrature;
  typedef detran_utilities::vec_int                     vec_int;
  typedef detran_utilities::vec_dbl                     vec_dbl;
  typedef detran_utilities::vec2_dbl                    vec2_dbl;
  typedef detran_utilities::vec3_dbl                    vec3_dbl;
  typedef detran_utilities::vec4_dbl                    vec4_dbl;
  typedef detran_utilities::vec5_dbl                    vec5_dbl;
  typedef detran_utilities::vec_size_t                  vec_size_t;
  typedef detran_utilities::vec2_size_t                 vec2_size_t;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /// Constructor
  ResponseSourceDetran(SP_node node, SP_indexer indexer);

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE
  //--------------------------------------------------------------------------//

  void compute(SP_response response, const ResponseIndex &index);

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  SP_solver solver() {return d_solver;}

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Detran db
  SP_db d_db;
  /// Detran material
  SP_material d_material;
  /// Detran mesh
  SP_mesh d_mesh;
  /// Detran solver
  SP_solver d_solver;
  /// Detran boundary
  SP_boundary d_B;
  /// Detran quadrature
  SP_quadrature d_quadrature;
  /// Energy basis [number surfaces]
  vec_basis d_basis_e;
  /// Spatial basis [number surfaces, dim per surface]
  vec2_basis d_basis_s;
  /// Azimuthal basis [number surfaces]
  vec_basis d_basis_a;
  /// Polar basis [number surfaces]
  vec_basis d_basis_p;
  /// Flag for expanding in angular flux (or the current)
  bool d_expand_angular_flux;
  /// Spatial dimensions in play for given axis
  vec2_size_t d_spatial_dim;


  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  /**
   *  @brief Set the boundary condition for a given incident order
   *  @param boundary   Reference to boundary flux container
   *  @param index      The response index that defines the boundary condition
   */
  void set_boundary(const ResponseIndex &index);

  /**
   *  @brief Expand the flux and boundary responses
   *  @param boundary   Pointer to boundary flux container
   *  @param response   Pointer to response container
   *  @param index      The response index that defines the boundary condition
   */
  void expand(SP_response response, const ResponseIndex &index);

  /// Expand the flux-responses
  void expand_boundary(SP_response response, const ResponseIndex &index);

  /// Construct the basis
  void construct_basis();

  /// Specialized basis construction routines
  //@{
  void construct_energy_basis();
  void construct_angular_basis_1D();
  void construct_angular_basis_2D();
  void construct_angular_basis_3D();
  //@}


  /// Sign to avoid integrating space variables in reverse.
  int spatial_sign(const size_t surface)
  {
    double sign = 1;
    if (surface == d_mesh->EAST  or
        surface == d_mesh->SOUTH or
        surface == d_mesh->BOTTOM)
    {
      sign = -1;
    }
    return sign;
  }

};

} // end namespace erme_response

#endif // erme_response_RESPONSESOURCEDETRAN_HH_

//----------------------------------------------------------------------------//
//              end of file ResponseSourceDetran.hh
//----------------------------------------------------------------------------//
