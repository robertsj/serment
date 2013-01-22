//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   PostProcess.hh
 *  @author robertsj
 *  @date   Jan 22, 2013
 *  @brief  PostProcess class definition.
 */
//---------------------------------------------------------------------------//

#ifndef erme_utils_POSTPROCESS_HH_
#define erme_utils_POSTPROCESS_HH_

#include "erme_geometry/NodeList.hh"
#include "erme_response/ResponseIndexer.hh"
#include "erme_response/ResponseServer.hh"
#include "erme/StateERME.hh"
#include "erme/ResponseMatrix.hh"
#include "erme/Connect.hh"
#include "erme/FissionOperator.hh"
#include "erme/AbsorptionOperator.hh"
#include "erme/LeakageOperator.hh"
#include "utilities/DBC.hh"
#include "utilities/SP.hh"
#include "utilities/Definitions.hh"

namespace erme_utils
{

/**
 *  @class PostProcess
 *  @brief Basic methods for processing converged results.
 */
class PostProcess
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<PostProcess>           SP_postprocess;
  typedef detran_utilities::InputDB::SP_input         SP_db;
  typedef erme::StateERME::SP_state                   SP_state;
  typedef erme::ResponseMatrix::SP_responsematrix     SP_R;
  typedef erme::Connect::SP_connect                   SP_M;
  typedef erme::FissionOperator::SP_fission           SP_F;
  typedef erme::AbsorptionOperator::SP_absorption     SP_A;
  typedef erme::LeakageOperator::SP_leakage           SP_L;
  typedef erme_geometry::NodeList::SP_nodelist        SP_nodelist;
  typedef erme_response::ResponseIndexer::SP_indexer  SP_indexer;
  typedef erme_response::ResponseServer::SP_server    SP_server;
  typedef detran_utilities::vec_dbl                   vec_dbl;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /// Constructor
  PostProcess(SP_db db,
              SP_nodelist nodes,
              SP_indexer indexer,
              SP_server server,
              SP_state state,
              SP_R R,
              SP_M M,
              SP_F F,
              SP_A A,
              SP_L L);

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  vec_dbl nodal_fission_density(const double norm = 1.0);

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Parameter database
  SP_db d_db;
  /// Node list
  SP_nodelist d_nodes;
  /// Indexer
  SP_indexer d_indexer;
  /// Response server
  SP_server d_server;
  /// State vector
  SP_state d_state;
  /// Response matrix
  SP_R d_R;
  /// Connectivity matrix
  SP_M d_M;
  /// Fission operator
  SP_F d_F;
  /// Absorption operator
  SP_A d_A;
  /// Leakage operator
  SP_L d_L;

};

} // end namespace erme_utils

#endif /* erme_utils_POSTPROCESS_HH_ */
