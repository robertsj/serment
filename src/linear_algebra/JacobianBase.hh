//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  JacobianBase.hh
 *  @brief JacobianBase class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//


#ifndef linear_algebra_JACOBIANBASE_HH_
#define linear_algebra_JACOBIANBASE_HH_

#include "MatrixBase.hh"
#include "utilities/SP.hh"

namespace linear_algebra
{

/**
 *  @class Jacobian
 *  @brief Defines the interface for Jacobian matrices
 */
class JacobianBase
{
public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<JacobianBase>    SP_jacobian;
  typedef Vector::SP_vector                     SP_vector;
  typedef MatrixBase::SP_matrix                 SP_matrix;

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// virtual destructor
  virtual ~JacobianBase(){}

  /// Update the Jacobian with the latest unknowns
  virtual void update(SP_vector x) = 0;

  virtual int local_size() const
  {
    return d_matrix->number_local_rows();
  }

  virtual int global_size() const
  {
    return d_matrix->number_global_rows();

  }

  SP_matrix matrix()
  {
    return d_matrix;
  }

protected:

  SP_matrix d_matrix;

};

} // end namespace linear_algebra

#endif /* linear_algebra_JACOBIANBASE_HH_ */

//----------------------------------------------------------------------------//
//                 end of JacobianBase.hh
//----------------------------------------------------------------------------//
