//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   AbsorptionOperator.hh
 * \brief  AbsorptionOperator 
 * \author Jeremy Roberts
 * \date   Aug 24, 2012
 */
//---------------------------------------------------------------------------//

#ifndef ABSORPTIONOPERATOR_HH_
#define ABSORPTIONOPERATOR_HH_

#include "ResponseOperator.hh"
#include "SP.hh"

namespace erme
{

class AbsorptionOperator: public ResponseOperator
{

public:

  typedef detran::SP<AbsorptionOperator> SP_absorption;

  AbsorptionOperator();

private:


};

} // end namespace detran

#endif // ABSORPTIONOPERATOR_HH_ 

//---------------------------------------------------------------------------//
//              end of file AbsorptionOperator.hh
//---------------------------------------------------------------------------//
