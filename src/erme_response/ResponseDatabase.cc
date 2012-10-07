//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ResponseDatabase.cc
 *  @author robertsj
 *  @date   Oct 1, 2012
 *  @brief  ResponseDatabase member definitions
 */
//---------------------------------------------------------------------------//

#include "ResponseDatabase.hh"

namespace erme_response
{

// Global instance
ResponseDatabase::SP_rfdb
ResponseDatabase::d_instance = ResponseDatabase::SP_rfdb(NULL);

//---------------------------------------------------------------------------//
ResponseDatabase::ResponseDatabase(std::string filename)
  : d_filename(filename)
  , d_open(false)
{

  // Open the HDF5 file
  d_file_id = H5Fopen(d_filename.c_str(),   // filename
                      H5F_ACC_RDONLY,       // read only
                      H5P_DEFAULT);         // file access property list
  Insist(d_file_id >= 0, "Error opening HDF5 file " + d_filename);
  d_open = true;

  //-------------------------------------------------------------------------//
  // LOAD DATABASE
  //-------------------------------------------------------------------------//

  // Get the number of nodes.  Each node is a group in the root group "/"
  herr_t ierr;
  ierr = H5Gget_num_objs(d_file_id, &d_number_nodes);
  Assert(!ierr);
  Insist(d_number_nodes, "The HDF5 file " + d_filename + " has no nodes!");
  std::cout << " # nodes = " << d_number_nodes << std::endl;

  // Loop over all nodes
  for (int i = 0; i < d_number_nodes; ++i)
  {
    // Get size of node name, add 1 for null terminator.
    ssize_t size = 1 + H5Lget_name_by_idx(d_file_id, ".",
                   H5_INDEX_NAME, H5_ITER_INC, i, NULL, 0, H5P_DEFAULT);

    // Allocate storage for name.
    char * name(new char[size]);

    // Retrieve name
    size = H5Lget_name_by_idx (d_file_id, ".",
           H5_INDEX_NAME, H5_ITER_INC, i, name, (size_t) size, H5P_DEFAULT);
    std::string nodename(name);
    delete [] name;
    std::cout << " READING NODE " << nodename << std::endl;

    // Open group
    hid_t group = H5Gopen(d_file_id, nodename.c_str(), H5P_DEFAULT);

    // Get the number of keff terms
    int number_keffs = 0;
    read_scalar_attribute(group, "number_keffs", number_keffs);
    std::cout << "   number_keffs " << number_keffs << std::endl;

    // Response size
    int response_size = 0;
    read_scalar_attribute(group, "response_size", response_size);
    std::cout << "   response_size " << response_size << std::endl;

    // Number of surfaces
    int number_surfaces = 0;
    read_scalar_attribute(group, "number_surfaces", number_surfaces);
    std::cout << "   number_surfaces " << number_surfaces << std::endl;

    // Number of surfaces
    int scheme = 0;
    read_scalar_attribute(group, "scheme", scheme);
    std::cout << "   scheme " << scheme << std::endl;

    // Instantiate this node.  This *assumes*  the db doesn't have repeats
    d_responses[nodename].number_keffs = number_keffs;
    d_responses[nodename].scheme = scheme;

    // Get keffs if needed
    if (scheme)
    {
      d_responses[nodename].keffs.resize(number_keffs, 0.0);
      hid_t dset = H5Dopen(group, "keffs", H5P_DEFAULT);
      herr_t status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                              H5P_DEFAULT, &d_responses[nodename].keffs[0]);
      status = H5Dclose(dset);

      for (int k = 0; k < number_keffs; ++k)
        std::cout << " k[" << k << "] = " << d_responses[nodename].keffs[k] << std::endl;
    }

    // Load the responses
    for (int k = 0; k < number_keffs; ++k)
    {
      SP_response r(new NodeResponse(response_size, number_surfaces));

      // Boundary response
      {
        hid_t dset  = H5Dopen(group, "R", H5P_DEFAULT);
        hid_t space = H5Dget_space(dset);
        hsize_t       dims_out[3];
        int rank    = H5Sget_simple_extent_ndims(space);
        int status  = H5Sget_simple_extent_dims(space, dims_out, NULL);

        // Want hyperslab corresponding to this k, i.e. data[k,:,:]
        hsize_t count[]  = {1, response_size, response_size};
        hsize_t offset[] = {k, 0, 0};
        status = H5Sselect_hyperslab(space, H5S_SELECT_SET, offset, NULL, count, NULL);

        // Define memory dataspace
        hsize_t dims[] = {response_size, response_size};
        double *buffer(new double[response_size*response_size]);
        hid_t memspace = H5Screate_simple (2, dims, NULL);

        // Read the set and kill the buffer
        status = H5Dread(dset, H5T_NATIVE_DOUBLE, memspace, space, H5P_DEFAULT, buffer);
        for (int ii = 0; ii < response_size; ++ii)
          for (int jj = 0; jj < response_size; ++jj)
            r->boundary_response(ii, jj) = buffer[ii + jj * response_size];
        delete [] buffer;

        H5Dclose(dset);
        H5Sclose(space);
      }

      // Leakage response
      {
        hid_t dset  = H5Dopen(group, "L", H5P_DEFAULT);
        hid_t space = H5Dget_space(dset);
        hsize_t       dims_out[3];
        int rank    = H5Sget_simple_extent_ndims(space);
        int status  = H5Sget_simple_extent_dims(space, dims_out, NULL);

        // Want hyperslab corresponding to this k, i.e. data[k,:,:]
        hsize_t count[]  = {1, number_surfaces, response_size};
        hsize_t offset[] = {k, 0, 0};
        status = H5Sselect_hyperslab(space, H5S_SELECT_SET, offset, NULL, count, NULL);

        // Define memory dataspace
        hsize_t dims[] = {number_surfaces, response_size};
        double *buffer(new double[number_surfaces*response_size]);
        hid_t memspace = H5Screate_simple (2, dims, NULL);

        // Read the set and kill the buffer
        status = H5Dread(dset, H5T_NATIVE_DOUBLE, memspace, space, H5P_DEFAULT, buffer);
        for (int ii = 0; ii < number_surfaces; ++ii)
          for (int jj = 0; jj < response_size; ++jj)
            r->leakage_response(ii, jj) = buffer[ii + jj * number_surfaces];
        delete [] buffer;

        H5Dclose(dset);
        H5Sclose(space);
      }

      // Fission response
      {
        hid_t dset  = H5Dopen(group, "F", H5P_DEFAULT);
        hid_t space = H5Dget_space(dset);
        hsize_t       dims_out[2];
        int rank    = H5Sget_simple_extent_ndims(space);
        int status  = H5Sget_simple_extent_dims(space, dims_out, NULL);

        // Want hyperslab corresponding to this k, i.e. data[k,:]
        hsize_t count[]  = {1, response_size};
        hsize_t offset[] = {k, 0};
        status = H5Sselect_hyperslab(space, H5S_SELECT_SET, offset, NULL, count, NULL);

        // Define memory dataspace
        hsize_t dims[] = {number_surfaces, response_size};
        hid_t memspace = H5Screate_simple (1, dims, NULL);

        // Read the set and kill the buffer
        status = H5Dread(dset, H5T_NATIVE_DOUBLE, memspace, space, H5P_DEFAULT,
                         &r->fission_response(0));

        H5Dclose(dset);
        H5Sclose(space);
      }

      // Absorption response
      {
        hid_t dset  = H5Dopen(group, "A", H5P_DEFAULT);
        hid_t space = H5Dget_space(dset);
        hsize_t       dims_out[2];
        int rank    = H5Sget_simple_extent_ndims(space);
        int status  = H5Sget_simple_extent_dims(space, dims_out, NULL);

        // Want hyperslab corresponding to this k, i.e. data[k,:]
        hsize_t count[]  = {1, response_size};
        hsize_t offset[] = {k, 0};
        status = H5Sselect_hyperslab(space, H5S_SELECT_SET, offset, NULL, count, NULL);

        // Define memory dataspace
        hsize_t dims[] = {number_surfaces, response_size};
        hid_t memspace = H5Screate_simple (1, dims, NULL);

        // Read the set and kill the buffer
        status = H5Dread(dset, H5T_NATIVE_DOUBLE, memspace, space, H5P_DEFAULT,
                         &r->absorption_response(0));

        H5Dclose(dset);
        H5Sclose(space);
      }

      // Add the temporary response to the vector
      d_responses[nodename].responses.push_back(r);

    } // end keff loop

    // Close group
    ierr = H5Gclose(group);

  } // end node loop

}

//---------------------------------------------------------------------------//
ResponseDatabase::~ResponseDatabase()
{
  // Close the HDF5 file \todo why does this yield an error? h5py related?
  //herr_t ierr = H5Fclose(d_file_id);
  //Ensure(!ierr);
  d_open = false;
}


//---------------------------------------------------------------------------//
ResponseDatabase::SP_rfdb ResponseDatabase::Create(std::string filename)
{
  //Preconditions
  Require(filename != "");

  // Create the instance if not done already
  if (!d_instance) d_instance = new ResponseDatabase(filename);

  // Postconditions
  Ensure(d_instance);
  return d_instance;
}

//---------------------------------------------------------------------------//
ResponseDatabase::SP_rfdb ResponseDatabase::Instance()
{
  Insist(d_instance, "Cannot use the database without it being created.")
  return d_instance;
}

} // end namespace erme_response


