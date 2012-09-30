//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   NodeResponse.hh
 * \brief  NodeResponse class definition.
 * \author Jeremy Roberts
 * \date   Aug 28, 2012
 */
//---------------------------------------------------------------------------//

#ifndef NODERESPONSE_HH_
#define NODERESPONSE_HH_

#include "utilities/DBC.hh"
#include "utilities/SP.hh"
#include <vector>

namespace erme_response
{

/*!
 *  \class NodeResponse
 *  \brief Container for nodal response functions
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
/*!
 *  \example erme_response/test/test_NodeResponse.cc
 *
 *  Test of NodeResponse class.
 */
class NodeResponse
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<NodeResponse>    SP_response;
  typedef unsigned int                          size_t;
  typedef std::vector<size_t>                   vec_size_t;
  typedef std::vector<double>                   vec_dbl;
  typedef std::vector<vec_dbl>                  vec2_dbl;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   *  \param moments_size   Size of node moments vector for each surface
   */
  NodeResponse(const size_t N, const size_t number_surfaces);

  //-------------------------------------------------------------------------//
  // ACCESS
  //-------------------------------------------------------------------------//

  /// Const access to boundary response
  const double& boundary_response(const size_t out,
                                  const size_t in) const;

  /// Mutable access to boundary response
  double& boundary_response(const size_t out,
                            const size_t in);

  /// Const access to fission response
  const double& fission_response(const size_t in) const;

  /// Mutable access to fission response
  double& fission_response(const size_t in);

  /// Const access to absorption response
  const double& absorption_response(const size_t in) const;

  /// Mutable access to absorption response
  double& absorption_response(const size_t in);

  /// Const access to leakage response
  const double& leakage_response(const size_t surface,
                                 const size_t in) const;

  /// Mutable access to leakage response
  double& leakage_response(const size_t surface,
                           const size_t in);

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

  /// Display the response data
  void display() const;

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Moment size
  const size_t d_N;

  /// Number of surfaces
  const size_t d_number_surfaces;

  /*!
   *  \brief Boundary function moments [N][N]
   *
   *  Stored [incoming][outgoing] so that for an incident condition,
   *  all outgoing values are contiguous in memory.
   */
  vec2_dbl d_boundary_response;

  /// Fission response [N]
  vec_dbl d_fission_response;

  /// Absorption response [N]
  vec_dbl d_absorption_response;

  /// Leakage response [nsurface][N]
  vec2_dbl d_leakage_response;

  /// Pin fission responses [Npins][N]
  vec2_dbl d_pin_response;

};

} // end namespace erme_response

//---------------------------------------------------------------------------//
// INLINE MEMBER DEFINITIONS
//---------------------------------------------------------------------------//

#include "NodeResponse.i.hh"

#endif // NODERESPONSE_HH_ 

//---------------------------------------------------------------------------//
//              end of file NodeResponse.hh
//---------------------------------------------------------------------------//
