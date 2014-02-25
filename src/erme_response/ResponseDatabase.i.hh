//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ResponseDatabase.i.hh
 *  @brief ResponseDatabase inline member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef erme_response_RESPONSEDATABASE_I_HH_
#define erme_response_RESPONSEDATABASE_I_HH_

#include "Interpolation.hh"

namespace erme_response
{

//----------------------------------------------------------------------------//
// \todo need to consider other interpolation and pin responses
inline void ResponseDatabase::get(std::string     nodename,
                                  SP_response     response,
                                  ResponseIndex   index,
                                  const double    keff)
{
  // Preconditions
  Require(nodename != "");
  Require(response);

  // Ensure we have the response.
  response_it it;
  it = d_responses.find(nodename);
  Insist(it != d_responses.end(), "node " + nodename
         + " not found in database " + d_filename);
  DBResponse &rf = it->second;

  // We are interpolating
  if (rf.scheme == 1)
  {
    // Determine number of points and size the vectors.
    size_t nk = rf.number_keffs;

    // If nk is larger than one, we need the bounding keffs.  For quadratic,
    // we want the central abscissa to be as close to the in keff as possible.
    if (nk > 1)
    {
      if (d_interpolation_order < nk + 1) nk = d_interpolation_order + 1;
    }
    vec_int kidx(nk, 0);
    vec_dbl kval(nk, 0.0);
    vec_dbl rval(nk, 0.0);

    if (nk > 1)
    {
      int k1 = nk / 2;
      for (; k1 < rf.number_keffs; ++k1)
        if (keff < rf.keffs[k1]) break;
      for (int i = 0; i < nk; ++i)
        kidx[i] = (i - nk/2) + k1;
      if (kidx[nk-1] >= rf.number_keffs)
      {
        for (int i = 0; i < nk; ++i)
          kidx[i] -= (kidx[nk-1] - rf.number_keffs + 1);
      }
    }
    // interpolation abscissa
    for (int i = 0; i < nk; ++i)
      kval[i] = rf.keffs[kidx[i]];

//    std::cout << " all keffs = " << std::endl;
//    for (int i = 0; i < rf.number_keffs; ++i)
//    {
//      std::cout << " k(i) = " << rf.keffs[i] << std::endl;
//    }
//    std::cout << " requested k = " << keff << std::endl;
//    std::cout << " bounding keffs = " << std::endl;
//    for (int i = 0; i < nk; ++i)
//    {
//      std::cout << " idx = " << kidx[i] << " k = " << kval[i] << std::endl;
//    }

    // incident nodal response index
    int in = index.nodal;

    // fill the responses
    for (int o = 0; o < response->size(); ++o)
    {
      Assertv(response->size() == rf.responses[kidx[0]]->size(),
              AsString(response->size()) + " vs " +
              AsString(rf.responses[kidx[0]]->size()));
      for (int i = 0; i < nk; ++i)
        rval[i] = rf.responses[kidx[i]]->boundary_response(o, in);
      response->boundary_response(o, in) = interpolate(keff, kval, rval);
    }
    for (int o = 0; o < response->number_surfaces(); ++o)
    {
      for (int i = 0; i < nk; ++i)
        rval[i] = rf.responses[kidx[i]]->leakage_response(o, in);
      response->leakage_response(o, in) = interpolate(keff, kval, rval);
    }
    for (int i = 0; i < nk; ++i)
      rval[i] = rf.responses[kidx[i]]->fission_response(in);
    response->fission_response(in) = interpolate(keff, kval, rval);
    for (int i = 0; i < nk; ++i)
      rval[i] = rf.responses[kidx[i]]->absorption_response(in);
    response->absorption_response(in) = interpolate(keff, kval, rval);
  }
  else
  {
    THROW("R(k) = R0 + R1/k ... EXPANSION NOT YET IMPLEMENTED");
  }

}

//----------------------------------------------------------------------------//
template <class T>
inline bool ResponseDatabase::read_scalar_attribute
(hid_t group, const char* name, T &value)
{
  detran_ioutils::HDF5_MemoryType mem;
  if (H5Aexists(group, name) == 0) return false;
  hid_t att = H5Aopen_name(group, name);
  herr_t status = H5Aread(att, mem.type<T>(), &value);
  Assert(!status);
  status = H5Aclose(att);
  Assert(!status);
  return true;
}

} // end namespace erme_response

#endif /* erme_response_RESPONSEDATABASE_I_HH_ */

//----------------------------------------------------------------------------//
//              end of file ResponseDatabase.i.hh
//----------------------------------------------------------------------------//
