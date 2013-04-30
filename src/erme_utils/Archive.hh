//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  Archive.hh
 *  @brief Archive class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef erme_utils_ARCHIVE_HH_
#define erme_utils_ARCHIVE_HH_

#include "utilities/InputDB.hh"
#include "erme_geometry/NodeList.hh"
#include <string>

namespace erme_utils
{

/**
 *  @class Archive
 *  @brief Use binary archives for reading and writing problem definition
 */
/**
 *  @example erme_utils/test/test_Archive.cc
 *  @brief   Test of Archive
 */
class Archive
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::InputDB::SP_input   SP_db;
  typedef erme_geometry::NodeList::SP_nodelist  SP_nodelist;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  Archive();

  virtual ~Archive();

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /// Return the database
  SP_db get_db();

  /// Return the node list
  SP_nodelist get_nodes();

  /**
   *  @brief Archive an input database and nodes
   *  @param db     User parameter database
   *  @param nodes  Complete list of nodes, not necessarily partitioned
   *  @param name   Name of the archive file
   */
  void archive(SP_db             db,
               SP_nodelist       nodes,
               const std::string filename = "serment.archive");

  /**
   *  @brief Unpack an archive
   *  @param name   Name of the archive file
   */
  void unarchive(const std::string filename);

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Parameter database
  SP_db d_db;
  /// Node list
  SP_nodelist d_nodes;

};

} /* namespace erme_utils */

#endif /* erme_utils_ARCHIVE_HH_ */
