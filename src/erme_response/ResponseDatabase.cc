//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ResponseDatabase.cc
 *  @brief ResponseDatabase member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "ResponseDatabase.hh"

namespace erme_response
{

// Global instance
ResponseDatabase::SP_rfdb
ResponseDatabase::d_instance = ResponseDatabase::SP_rfdb(NULL);

//----------------------------------------------------------------------------//
ResponseDatabase::ResponseDatabase(std::string filename, size_t order, std::string type)
  : d_filename(filename)
  , d_open(false)
  , d_interpolation_order(order)
  , d_interpolation_type(type)
{
  bool db = false;

  // Open the HDF5 file
  d_file_id = H5Fopen(d_filename.c_str(),   // filename
                      H5F_ACC_RDONLY,       // read only
                      H5P_DEFAULT);         // file access property list
  Insist(d_file_id >= 0, "Error opening HDF5 file " + d_filename);
  d_open = true;

  //--------------------------------------------------------------------------//
  // LOAD DATABASE
  //--------------------------------------------------------------------------//

  // Get the number of nodes.  Each node is a group in the root group "/"
  herr_t ierr;
  ierr = H5Gget_num_objs(d_file_id, &d_number_nodes);
  Assert(!ierr);
  Insist(d_number_nodes, "The HDF5 file " + d_filename + " has no nodes!");

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
    if (db) std::cout << " READING NODE " << nodename << std::endl;

    // Open group
    hid_t group = H5Gopen(d_file_id, nodename.c_str(), H5P_DEFAULT);

    // Read attributes
    bool flag = false;

    // Get the number of keff terms
    int number_keffs = 0;
    flag = read_scalar_attribute(group, "number_keffs", number_keffs);
    if (db) std::cout << "   number_keffs " << number_keffs << std::endl;

    // Response size
    int response_size = 0;
    flag = read_scalar_attribute(group, "response_size", response_size);
    if (db) std::cout << "   response_size " << response_size << std::endl;

    // Number of surfaces
    int number_surfaces = 0;
    flag = read_scalar_attribute(group, "number_surfaces", number_surfaces);
    if (db) std::cout << "   number_surfaces " << number_surfaces << std::endl;

    // Number of pins
    int number_pins = 0;
    flag = read_scalar_attribute(group, "number_pins", number_pins);
    if (db) std::cout << "   number_pins " << number_pins << std::endl;

    // Response expansion scheme
    int scheme = 0;
    flag = read_scalar_attribute(group, "scheme", scheme);
    if (db) std::cout << "   scheme " << scheme << std::endl;

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
      Assert(!status);

      if (db)
      {
        for (int k = 0; k < number_keffs; ++k)
          std::cout << " k[" << k << "] = " << d_responses[nodename].keffs[k] << std::endl;
      }
    }

    // Load the responses
    for (int k = 0; k < number_keffs; ++k)
    {
      SP_response r(new NodeResponse(response_size, number_surfaces,
                                     number_pins));

      // Boundary response
      {
        hid_t dset  = H5Dopen(group, "R", H5P_DEFAULT);
        hid_t space = H5Dget_space(dset);
        hsize_t       dims_out[3];
        int rank    = H5Sget_simple_extent_ndims(space);
        int ndims   = H5Sget_simple_extent_dims(space, dims_out, NULL);
        Assert(ndims > 0);

        // Want hyperslab corresponding to this k, i.e. data[k,:,:]
        hsize_t count[]  = {1, response_size, response_size};
        hsize_t offset[] = {k, 0, 0};
        int status = H5Sselect_hyperslab(space, H5S_SELECT_SET, offset, NULL, count, NULL);
        Assert(!status);

        // Define memory dataspace
        hsize_t dims[] = {response_size, response_size};
        double *buffer(new double[response_size*response_size]);
        hid_t memspace = H5Screate_simple (2, dims, NULL);

        // Read the set and kill the buffer
        status = H5Dread(dset, H5T_NATIVE_DOUBLE, memspace, space, H5P_DEFAULT, buffer);
        Assert(!status);
        // Note the indexing of ii and jj: jj changes fastest.
        for (int ii = 0; ii < response_size; ++ii)
          for (int jj = 0; jj < response_size; ++jj)
            r->boundary_response(ii, jj) = buffer[jj + ii * response_size];
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
        int ndims   = H5Sget_simple_extent_dims(space, dims_out, NULL);
        Assert(ndims > 0);

        // Want hyperslab corresponding to this k, i.e. data[k,:,:]
        hsize_t count[]  = {1, number_surfaces, response_size};
        hsize_t offset[] = {k, 0, 0};
        int status = H5Sselect_hyperslab(space, H5S_SELECT_SET, offset, NULL, count, NULL);
        Assert(!status);

        // Define memory dataspace
        hsize_t dims[] = {number_surfaces, response_size};
        double *buffer(new double[number_surfaces*response_size]);
        hid_t memspace = H5Screate_simple (2, dims, NULL);

        // Read the set and kill the buffer
        status = H5Dread(dset, H5T_NATIVE_DOUBLE, memspace, space, H5P_DEFAULT, buffer);
        Assert(!status);
//        for (int kk = 0; kk < number_surfaces*response_size; ++kk)
//        {
//          std::cout << " buffer " << kk << " = " << buffer[kk] << std::endl;
//        }
        for (int ii = 0; ii < number_surfaces; ++ii)
          for (int jj = 0; jj < response_size; ++jj)
            r->leakage_response(ii, jj) = buffer[jj + ii * response_size];
        delete [] buffer;

        H5Dclose(dset);
        H5Sclose(space);
      }

      // Pin powers
      if (number_pins > 0)
      {
        hid_t dset  = H5Dopen(group, "pin_power", H5P_DEFAULT);
        hid_t space = H5Dget_space(dset);
        hsize_t       dims_out[3];
        int rank    = H5Sget_simple_extent_ndims(space);
        int ndims   = H5Sget_simple_extent_dims(space, dims_out, NULL);
        Assert(ndims > 0);

        // Want hyperslab corresponding to this k, i.e. data[k,:,:]
        hsize_t count[]  = {1, number_pins, response_size};
        hsize_t offset[] = {k, 0, 0};
        int status = H5Sselect_hyperslab(space, H5S_SELECT_SET, offset, NULL, count, NULL);
        Assert(!status);

        // Define memory dataspace
        hsize_t dims[] = {number_pins, response_size};
        double *buffer(new double[number_pins*response_size]);
        hid_t memspace = H5Screate_simple (2, dims, NULL);

        // Read the set and kill the buffer
        status = H5Dread(dset, H5T_NATIVE_DOUBLE, memspace, space, H5P_DEFAULT, buffer);
        Assert(!status);

        for (int ii = 0; ii < number_pins; ++ii)
          for (int jj = 0; jj < response_size; ++jj)
            r->pin_power(ii, jj) = buffer[jj + ii * response_size];
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
        int ndims   = H5Sget_simple_extent_dims(space, dims_out, NULL);
        Assert(ndims > 0);

        // Want hyperslab corresponding to this k, i.e. data[k,:]
        hsize_t count[]  = {1, response_size};
        hsize_t offset[] = {k, 0};
        int status = H5Sselect_hyperslab(space, H5S_SELECT_SET, offset, NULL, count, NULL);
        Assert(!status);

        // Define memory dataspace
        hsize_t dims[] = {response_size};
        hid_t memspace = H5Screate_simple (1, dims, NULL);

        // Read the set and kill the buffer
        status = H5Dread(dset, H5T_NATIVE_DOUBLE, memspace, space, H5P_DEFAULT,
                         &r->fission_response(0));
        Assert(!status);

        H5Dclose(dset);
        H5Sclose(space);
      }

      // Absorption response
      {
        hid_t dset  = H5Dopen(group, "A", H5P_DEFAULT);
        hid_t space = H5Dget_space(dset);
        hsize_t       dims_out[2];
        int rank    = H5Sget_simple_extent_ndims(space);
        int ndims   = H5Sget_simple_extent_dims(space, dims_out, NULL);
        Assert(ndims > 0);

        // Want hyperslab corresponding to this k, i.e. data[k,:]
        hsize_t count[]  = {1, response_size};
        hsize_t offset[] = {k, 0};
        int status = H5Sselect_hyperslab(space, H5S_SELECT_SET, offset, NULL, count, NULL);
        Assert(!status);

        // Define memory dataspace
        hsize_t dims[] = {response_size};
        hid_t memspace = H5Screate_simple (1, dims, NULL);

        // Read the set and kill the buffer
        status = H5Dread(dset, H5T_NATIVE_DOUBLE, memspace, space, H5P_DEFAULT,
                         &r->absorption_response(0));
        Assert(!status);

        H5Dclose(dset);
        H5Sclose(space);
      }

      // Nodal power
      {
        hid_t dset  = H5Dopen(group, "nodal_power", H5P_DEFAULT);
        hid_t space = H5Dget_space(dset);
        hsize_t       dims_out[2];
        int rank    = H5Sget_simple_extent_ndims(space);
        int ndims   = H5Sget_simple_extent_dims(space, dims_out, NULL);
        Assert(ndims > 0);

        // Want hyperslab corresponding to this k, i.e. data[k,:]
        hsize_t count[]  = {1, response_size};
        hsize_t offset[] = {k, 0};
        int status = H5Sselect_hyperslab(space, H5S_SELECT_SET, offset, NULL, count, NULL);
        Assert(!status);

        // Define memory dataspace
        hsize_t dims[] = {response_size};
        hid_t memspace = H5Screate_simple (1, dims, NULL);

        // Read the set and kill the buffer
        status = H5Dread(dset, H5T_NATIVE_DOUBLE, memspace, space, H5P_DEFAULT,
                         &r->nodal_power(0));
        Assert(!status);

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

//----------------------------------------------------------------------------//
ResponseDatabase::~ResponseDatabase()
{
  // Close the HDF5 file \todo why does this yield an error? h5py related?
  //herr_t ierr = H5Fclose(d_file_id);
  //Ensure(!ierr);
  d_open = false;
}

//----------------------------------------------------------------------------//
ResponseDatabase::SP_rfdb
ResponseDatabase::Create(std::string filename, size_t order, std::string type)
{
  //Preconditions
  Require(filename != "");

  // Create the instance if not done already
  d_instance = new ResponseDatabase(filename, order, type);

  // Postconditions
  Ensure(d_instance);
  return d_instance;
}

//----------------------------------------------------------------------------//
ResponseDatabase::SP_rfdb ResponseDatabase::Instance()
{
  Insist(d_instance, "Cannot use the database without it being created.")
  return d_instance;
}

} // end namespace erme_response

//----------------------------------------------------------------------------//
//              end of file ResponseDatabase.cc
//----------------------------------------------------------------------------//

