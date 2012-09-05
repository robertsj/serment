//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   FissionOperator.cc
 * \author robertsj
 * \date   Sep 4, 2012
 * \brief  FissionOperator class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#include "FissionOperator.hh"

namespace erme
{

FissionOperator::FissionOperator(SP_nodelist nodes,
                                 SP_indexer indexer,
                                 SP_server server)
  : ResponseOperator(nodes, indexer, server)
  , Vector(indexer->number_local_moments())
{
  /* ... */
}

void FissionOperator::update()
{
  using std::cout;
  using std::endl;

  // Offset for a block.  Starts at this vector's lower bound.
  int offset = lower_bound();

  // Loop through nodes
  for (int n = 0; n < d_nodes->number_local_nodes(); n++)
  {

    // Get response
    SP_response r = d_server->response(n);

    // Row indices
    std::vector<int> indices(r->size(), offset);
    for (size_type i = 0; i < indices.size(); i++)
      indices[i] += i;

    // Insert the values (size, rows, values)
    insert_values(r->size(), &indices[0], &r->fission_response(0));

    // Increment the offset
    offset += indices.size();
  }

  // Assemble after finishing
  assemble();

}

} // end namespace erme

//---------------------------------------------------------------------------//
//              end of file FissionOperator.cc
//---------------------------------------------------------------------------//



