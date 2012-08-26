//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   FissionOperator.hh
 * \brief  FissionOperator 
 * \author Jeremy Roberts
 * \date   Aug 24, 2012
 */
//---------------------------------------------------------------------------//

#ifndef FISSIONOPERATOR_HH_
#define FISSIONOPERATOR_HH_

#include "ResponseOperator.hh"

namespace erme
{

class FissionOperator: public ResponseOperator
{

public:

  typedef detran::SP<FissionOperator> SP_fission;

  FissionOperator();

private:


};

} // end namespace detran

#endif // FISSIONOPERATOR_HH_ 

//---------------------------------------------------------------------------//
//              end of file FissionOperator.hh
//---------------------------------------------------------------------------//
