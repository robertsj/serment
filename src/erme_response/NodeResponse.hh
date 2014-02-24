//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  NodeResponse.hh
 *  @brief NodeResponse class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef erme_response_NODERESPONSE_HH_
#define erme_response_NODERESPONSE_HH_

#include "utilities/DBC.hh"
#include "utilities/SP.hh"
#include <vector>

namespace erme_response
{

/**
 *  @class NodeResponse
 *  @brief Container for nodal response functions
 *
 *  NodeResponse is a very simple container for responses.  For the simplest
 *  problems, the only responses required are for the boundary function
 *  (usually the partial current or angular flux), fission, absorption, and
 *  leakage.  Suppose we have a node with \f$ S \f$ surfaces, and on each
 *  surface, we have a partial current expanded in \f$ M \f$ terms.  The
 *  total boundary vector for this node has a size of \f$ N = SM \f$.  The
 *  total size of the boundary response data is then \f$ N^2 \f$, since
 *  each of the \f$ N \f$ moments contributes to outgoing currents that
 *  are also expanded in \f$ N \f$ terms.  The fission and absorption
 *  operators represent vectors with which the node incident boundary
 *  function is folded to yield total fission and absorption rates; that
 *  vector has a size of \f$ N \f$.  The leakage operator yields the
 *  total leakage from each of \f$ S \f$ when operated on the node
 *  boundary function; hence, the node leakage data has a size of
 *  \f$ SN \f$.
 *
 *  To be efficient, we want these responses stored contiguously if
 *  possible.  If we view the node boundary responses as
 *  an \f$ N \times N \f$ block, then one column is produced per
 *  incident response.  Hence, we are best served using column-oriented
 *  storage for this data.  Moreover, this facilitates the probable
 *  case of response server processes sending these columns back
 *  to a response driver process, after which the driver participates
 *  in the global solve.
 *
 *  In this initial implementation, we'll use STL vectors.  For the
 *  boundary responses, columns will be stored contiguously.
 *
 */
/**
 *  @example erme_response/test/test_NodeResponse.cc
 *
 *  Test of NodeResponse class.
 */
class NodeResponse
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<NodeResponse>    SP_response;
  typedef unsigned int                          size_t;
  typedef std::vector<size_t>                   vec_size_t;
  typedef std::vector<double>                   vec_dbl;
  typedef std::vector<vec_dbl>                  vec2_dbl;

  //--------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param N                Total number of nodal surface moments
   *  @param number_surfaces  Number of nodal surfaces
   *  @param number_pins      Number of pins in the node
   */
  NodeResponse(const size_t N,
               const size_t number_surfaces,
               const size_t number_pins = 0);

  //--------------------------------------------------------------------------//
  // ACCESS
  //--------------------------------------------------------------------------//

  /// Const access to boundary response
  double boundary_response(const size_t out, const size_t in) const;
  /// Mutable access to boundary response
  double& boundary_response(const size_t out, const size_t in);

  /// Const access to fission response
  double fission_response(const size_t in) const;
  /// Mutable access to fission response
  double& fission_response(const size_t in);

  /// Const access to absorption response
  double absorption_response(const size_t in) const;
  /// Mutable access to absorption response
  double& absorption_response(const size_t in);

  /// Const access to leakage response
  double leakage_response(const size_t surface, const size_t in) const;
  /// Mutable access to leakage response
  double& leakage_response(const size_t surface, const size_t in);

  /// Const access to leakage response
  double nodal_power(const size_t in) const;
  /// Mutable access to leakage response
  double& nodal_power(const size_t in);

  /// Const access to leakage response
  double pin_power(const size_t p, const size_t in) const;
  /// Mutable access to leakage response
  double& pin_power(const size_t p, const size_t in);

  /// Clear all the responses
  void clear();

  /// Return moment size
  size_t size() const
  {
    return d_N;
  }

  /// Return number of surfaces
  size_t number_surfaces() const
  {
    return d_number_surfaces;
  }

  /// Return number of pins
  size_t number_pins() const
  {
    return d_number_pins;
  }

  /// Display the response data
  void display() const;

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Moment size
  const size_t d_N;
  /// Number of surfaces
  const size_t d_number_surfaces;
  /// Number of pins
  const size_t d_number_pins;
  /// Boundary function moments [N][N]; outgoing is stored contiguously
  vec2_dbl d_boundary_response;
  /// Fission response [N]
  vec_dbl d_fission_response;
  /// Absorption response [N]
  vec_dbl d_absorption_response;
  /// Leakage response [nsurface][N]
  vec2_dbl d_leakage_response;
  /// Nodal power responses [N]
  vec_dbl d_nodal_power;
  /// Pin power response [Npins][N]
  vec2_dbl d_pin_power;

};

} // end namespace erme_response

//----------------------------------------------------------------------------//
// INLINE MEMBER DEFINITIONS
//----------------------------------------------------------------------------//

#include "NodeResponse.i.hh"

#endif // erme_response_NODERESPONSE_HH_

//----------------------------------------------------------------------------//
//              end of file NodeResponse.hh
//----------------------------------------------------------------------------//
