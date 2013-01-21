//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ResponseSourceDetran.hh
 *  @brief  ResponseSourceDetran class definition
 *  @author Jeremy Roberts
 *  @date   Sep 1, 2012
 */
//---------------------------------------------------------------------------//

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

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

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
  typedef detran_utilities::vec_size_t                  vec_size_t;
  typedef detran_utilities::vec2_size_t                 vec2_size_t;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /// Constructor
  ResponseSourceDetran(SP_node node, SP_indexer indexer);

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE
  //-------------------------------------------------------------------------//

  void compute(SP_response response, const ResponseIndex &index);

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  SP_solver solver() {return d_solver;}

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Detran db
  SP_db d_db;
  /// Detran material
  SP_material d_material;
  /// Detran mesh
  SP_mesh d_mesh;
  /// Detran solver
  SP_solver d_solver;
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
  /// Spatial dimensions in play for given axis
  vec2_size_t d_spatial_dim;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /**
   *  @brief Set the boundary condition for a given incident order
   *  @param boundary   Reference to boundary flux container
   *  @param index      The response index that defines the boundary condition
   */
  void set_boundary(B& boundary, const ResponseIndex &index);

  /**
   *  @brief Expand the flux and boundary responses
   *  @param boundary   Pointer to boundary flux container
   *  @param response   Pointer to response container
   *  @param index      The response index that defines the boundary condition
   */
  void expand(const B& boundary, SP_response response, const ResponseIndex &index);

  /// Expand the flux-responses
  void expand_flux(SP_response response, const ResponseIndex &index);

  /// Construct the basis
  void construct_basis();

};

} // end namespace erme_response

#endif // erme_response_RESPONSESOURCEDETRAN_HH_

//---------------------------------------------------------------------------//
//              end of file ResponseSourceDetran.hh
//---------------------------------------------------------------------------//
