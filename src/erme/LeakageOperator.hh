//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LeakageOperator.hh
 * \brief  LeakageOperator 
 * \author Jeremy Roberts
 * \date   Aug 24, 2012
 */
//---------------------------------------------------------------------------//

#ifndef LEAKAGEOPERATOR_HH_
#define LEAKAGEOPERATOR_HH_

#include "ResponseOperator.hh"

class LeakageOperator: public ResponseOperator
{

public:

  typedef detran::SP<LeakageOperator> SP_leakage;

  LeakageOperator();

private:


};

#endif // LEAKAGEOPERATOR_HH_ 

//---------------------------------------------------------------------------//
//              end of file LeakageOperator.hh
//---------------------------------------------------------------------------//
