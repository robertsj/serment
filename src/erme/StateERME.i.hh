//----------------------------------*-C++-*----------------------------------//
/**
 * @file   StateERME.i.hh
 * @brief  StateERME inline member definitions
 * @author Jeremy Roberts
 * @date   Aug 23, 2012
 */
//---------------------------------------------------------------------------//

#ifndef erme_STATEERME_I_HH_
#define erme_STATEERME_I_HH_

namespace erme
{

inline void StateERME::set_k(const double k_val)
{
  d_k = k_val;
}

inline void StateERME::set_lambda(const double lambda_val)
{
  d_lambda = lambda_val;
}

inline double StateERME::k() const
{
  return d_k;
}

inline double StateERME::lambda() const
{
  return d_lambda;
}

inline StateERME::size_t StateERME::local_size() const
{
  std::cout << "LOCALSIZE=" << d_local_size << std::endl;
  return d_local_size;
}

inline StateERME::size_t StateERME::global_size() const
{
  std::cout << "GLOBALSIZE=" << d_global_size << std::endl;

  return d_global_size;
}

inline const StateERME::Vector& StateERME::moments() const
{
  return d_boundary_moments;
}

inline StateERME::Vector& StateERME::moments()
{
  return d_boundary_moments;
}

} // end namespace erme

#endif // erme_STATEERME_I_HH_

//---------------------------------------------------------------------------//
//              end of file StateERME.i.hh
//---------------------------------------------------------------------------//
