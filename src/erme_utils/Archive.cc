//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Archive.cc
 *  @brief Archive member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "Archive.hh"
#include "comm/Comm.hh"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <fstream>

namespace erme_utils
{

//----------------------------------------------------------------------------//
Archive::Archive()
  : d_db(NULL)
{
  // TODO Auto-generated constructor stub
}

//----------------------------------------------------------------------------//
Archive::~Archive()
{
  // TODO Auto-generated destructor stub
}

//----------------------------------------------------------------------------//
Archive::SP_db Archive::get_db()
{
  if (serment_comm::Comm::world_rank() == 0)
    Require(d_db);
  return d_db;
}

//----------------------------------------------------------------------------//
Archive::SP_nodelist Archive::get_nodes()
{
  if (serment_comm::Comm::world_rank() == 0)
    Require(d_nodes);
  return d_nodes;
}

//----------------------------------------------------------------------------//
void Archive::archive(SP_db &db, SP_nodelist &nodes, const std::string filename)
{
  if (serment_comm::Comm::world_rank() == 0)
  {
    Insist(db,             "Database is NULL.");
    Insist(nodes,          "Node list is NULL.");
    Insist(filename != "", "Filename is empty.");

    std::ofstream outstream(filename.c_str());
    boost::archive::binary_oarchive outarchive(outstream);
    outarchive << db;
    outarchive << nodes;
    outstream.close();
  }
}

//----------------------------------------------------------------------------//
void Archive::unarchive(const std::string filename)
{
  if (serment_comm::Comm::world_rank() == 0)
  {
    Insist(filename != "", "Filename is empty.");

    std::ifstream instream(filename.c_str());
    boost::archive::binary_iarchive inarchive(instream);
    inarchive >> d_db;
    inarchive >> d_nodes;
    instream.close();

    Ensure(d_db);
    Ensure(d_nodes);
  }
}

void Archive::unarchive(const std::string filename, SP_db &db, SP_nodelist &nodes)
{
  if (serment_comm::Comm::world_rank() == 0)
  {
    Insist(filename != "", "Filename is empty.");

    std::ifstream instream(filename.c_str());
    boost::archive::binary_iarchive inarchive(instream);
    inarchive >> db;
    inarchive >> nodes;
    instream.close();

    Ensure(db);
    Ensure(nodes);
  }
}

} /* namespace erme_utils */

//----------------------------------------------------------------------------//
//              end of Archive.cc
//----------------------------------------------------------------------------//
