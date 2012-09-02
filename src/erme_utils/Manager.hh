//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Manager.hh
 * \brief  Manager class definition
 * \author Jeremy Roberts
 * \date   Sep 1, 2012
 */
//---------------------------------------------------------------------------//

#ifndef MANAGER_HH_
#define MANAGER_HH_

//---------------------------------------------------------------------------//
/*! \mainpage Serment: An Eigenvalue Response Matrix Code
 *
 * \section sec_introduction Introduction
 *
 * The eigenvalue response matrix method (ERMM) is a numerical approach
 * for static reactor analysis.  The basic idea behind the method is
 * to decompose a global domain in space and link the resulting independent
 * \e nodes via approximate boundary conditions.  Because the nodes are
 * completely independent in a computational sense, the "physical"
 * decomposition results in natural computational domain decomposition
 * for parallel computation.
 *
 * Since nodes communicate only at the boundaries,
 * a variety of transport methods, both deterministic and stochastic, can
 * be used to generate the boundary conditions as long as a consistent
 * approximation is used (<em>e.g.</em> a P1 approximation, etc.).
 * Furthermore, different approximations can be used for different nodes,
 * and varying levels of boundary approximation can be used throughout
 * the domain.
 *
 * A more detailed discussion of the methods used in Serment can be found
 * in the theory documentation.  Details on implementation specifics can
 * be found throughout the rest of this documentation.
 *
 *
 * \section install_sec Installation
 *
 * \subsection step1 Step 1: Opening the box
 *
 * etc...
 */
//---------------------------------------------------------------------------//

/*!
 *  \brief Namespace for higher level routines for organizing response
 *         matrix problems and their solutions
 */
namespace erme_utils
{


} // end namespace erme_utils

#endif // MANAGER_HH_ 

//---------------------------------------------------------------------------//
//              end of file Manager.hh
//---------------------------------------------------------------------------//
