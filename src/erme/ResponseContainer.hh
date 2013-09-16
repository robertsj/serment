//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ResponseContainer.hh
 *  @brief ResponseContainer class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef erme_RESPONSECONTAINER_HH_
#define erme_RESPONSECONTAINER_HH_

#include "erme/ResponseMatrix.hh"
#include "erme/Connect.hh"
#include "erme/FissionOperator.hh"
#include "erme/AbsorptionOperator.hh"
#include "erme/LeakageOperator.hh"

namespace erme
{

/**
 *  @struct ResponseContainer
 *  @brief  Convenience container for erme objects
 *
 *  This will make adding other responses (e.g. region power)
 *  less invasive.
 */
struct ResponseContainer
{

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

	typedef detran_utilities::SP<ResponseContainer>     SP_responsecontainer;
  typedef erme_geometry::NodeList::SP_nodelist        SP_nodelist;
  typedef erme_response::ResponseIndexer::SP_indexer  SP_indexer;
  typedef erme_response::ResponseServer::SP_server    SP_server;
  typedef ResponseMatrix::SP_responsematrix     			SP_R;
  typedef Connect::SP_connect                   			SP_M;
  typedef FissionOperator::SP_fission           			SP_F;
  typedef AbsorptionOperator::SP_absorption     			SP_A;
  typedef LeakageOperator::SP_leakage           			SP_L;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

	ResponseContainer(SP_nodelist nodes, SP_indexer indexer, SP_server server)
  {
		Insist(serment_comm::communicator == serment_comm::world,
				   "Responses should be built from the world communicator.");
		Require(nodes);
		Require(indexer);
		Require(server);

		if (serment_comm::Comm::is_global())
		{
		  std::cout << "********* BUILDING M ....";
			M = new Connect(nodes, indexer);
			std::cout << " done. " << std::endl;

			std::cout << "********* BUILDING R ....";
			R = new ResponseMatrix(nodes, indexer, server);
			std::cout << " done. " << std::endl;

			std::cout << "********* BUILDING L ....";
			L = new LeakageOperator(nodes, indexer, server);
			std::cout << " done. " << std::endl;

			std::cout << "********* BUILDING F ....";
			F = new FissionOperator(nodes, indexer, server);
			std::cout << " done. " << std::endl;

			std::cout << "********* BUILDING A ....";
			A = new AbsorptionOperator(nodes, indexer, server);
			std::cout << " done. " << std::endl;
		}
  }

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

	/// Update the response operators with the latest eigenvalue held by server
	void update()
	{
    if (serment_comm::Comm::is_global())
    {
      R->update();
      L->update();
      F->update();
      A->update();
    }
	}

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Response matrix
  SP_R R;
  /// Connectivity matrix
  SP_M M;
  /// Fission operator
  SP_F F;
  /// Absorption operator
  SP_A A;
  /// Leakage operator
  SP_L L;

};

} // end namespace erme


#endif /* erme_RESPONSECONTAINER_HH_ */
