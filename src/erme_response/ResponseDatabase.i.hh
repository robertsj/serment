//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ResponseDatabase.i.hh
 *  @author robertsj
 *  @date   Oct 1, 2012
 *  @brief  ResponseDatabase inline member definitions
 */
//---------------------------------------------------------------------------//

#ifndef RESPONSEDATABASE_I_HH_
#define RESPONSEDATABASE_I_HH_

namespace erme_response
{

// \todo need to consider other interpolation and pin responses
inline void ResponseDatabase::get(std::string     nodename,
                                  SP_response     response,
                                  ResponseIndex   index,
                                  const double    keff)
{
  // Preconditions
  Require(!nodename.compare(""));
  Require(response);

  // Ensure we have the response.
  response_it it;
  it = d_responses.find(nodename);
  Insist(it != d_responses.end(), "node " + nodename
         + " not found in datbase " + d_filename);
  DBResponse &rf = it->second;

  // We are interpolating
  //
  // Linear:    r(k) ~ [r2-r1]/[k2-k1] * k + [r2x1-r1x2]/[x1-x2]
  // Quadratic: r(k) ~ ... later
  //
  if (rf.scheme == 1)
  {
    int k0 = 0; // first k point
    int k1 = 0; // second k point
    if (rf.number_keffs > 1)
    {
      k1 = 1;
      for (; k1 < rf.number_keffs; ++k1)
        if (keff < rf.keffs[k1]) break;
      if (k1 > 1)  k0 = k1 - 1;
    }

    // incident nodal response index
    int in = index.nodal;
    // interpolation function values and abscissa
    double r0 = 0;
    double r1 = 0;
    double x0 = rf.keffs[k0];
    double x1 = rf.keffs[k1];
    // fill the responses
    for (int o = 0; o < response->size(); ++o)
    {
      r0 = rf.responses[k0]->boundary_response(o, in);
      r1 = rf.responses[k1]->boundary_response(o, in);
      response->boundary_response(o, in) = interpolate(keff, x0, x1, r0, r1);
    }
    for (int o = 0; o < response->number_surfaces(); ++o)
    {
      r0 = rf.responses[k0]->leakage_response(o, in);
      r1 = rf.responses[k1]->leakage_response(o, in);
      response->leakage_response(o, in) = interpolate(keff, x0, x1, r0, r1);
    }
    r0 = rf.responses[k0]->fission_response(in);
    r1 = rf.responses[k1]->fission_response(in);
    response->fission_response(in) = interpolate(keff, x0, x1, r0, r1);
    r0 = rf.responses[k0]->absorption_response(in);
    r1 = rf.responses[k1]->absorption_response(in);
    response->absorption_response(in) = interpolate(keff, x0, x1, r0, r1);
  }
  else
  {
    THROW("R(k) = R0 + R1/k ... EXPANSION NOT YET IMPLEMENTED");
  }

}

template <class T>
inline bool ResponseDatabase::read_scalar_attribute
(hid_t group, const char* name, T &value)
{
  detran_ioutils::HDF5_MemoryType mem;
  hid_t att     = H5Aopen_name(group, name);
  herr_t status = H5Aread(att, mem.type<T>(), &value);
  status        = H5Aclose(att);
  return true;
}

inline double ResponseDatabase::
interpolate(double x, double x0, double x1, double r0, double r1)
{
  double idx = 1.0 / (x1 - x0);
  return idx * ( (r1 - r0) * x + (r1 * x0 - r0 * x1) );
}

} // end namespace erme_response


#endif /* RESPONSEDATABASE_I_HH_ */
