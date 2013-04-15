//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Jacobian.cc
 *  @author robertsj
 *  @date   Apr 12, 2013
 *  @brief  Jacobian class definition.
 */
//---------------------------------------------------------------------------//

#include "Jacobian.hh"

namespace erme_solver
{

//---------------------------------------------------------------------------//
Jacobian::Jacobian(SP_state state,
									 SP_R 		R,
									 SP_M 		M,
									 SP_F 		F,
									 SP_A 		A,
									 SP_L 		L)
{

}

//---------------------------------------------------------------------------//
PetscErrorCode Jacobian::shell_multiply(Vec x, Vec y)
{
	PetscErrorCode ierr = 0;


	return ierr;
}

} // end namespace erme_solver

