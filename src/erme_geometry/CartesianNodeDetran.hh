//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CartesianNodeDetran.hh
 * \brief  CartesianNodeDetran class definition
 * \author Jeremy Roberts
 * \date   Aug 22, 2012
 */
//---------------------------------------------------------------------------//

#ifndef CARTESIANNODEDETRAN_HH_
#define CARTESIANNODEDETRAN_HH_

// ERME geometry
#include "CartesianNode.hh"

// Detran
#include "InputDB.hh"
#include "Material.hh"
#include "detran_geometry.hh"

namespace erme_geometry
{

/*!
 *  \class CartesianNodeDetran
 *  \brief Specialization of CartesianNode for use with Detran
 *
 *  The node essentially defines the local problem to be solved by
 *  Detran.  The required input consists of the parameters database,
 *  the material definitions, and the mesh.  Any other needed objects are
 *  created later in the response generator.
 *
 */
class CartesianNodeDetran: public CartesianNode
{

public:

  /// Useful typedefs
  typedef CartesianNode                   Base;
  typedef detran::InputDB::SP_input       SP_db;
  typedef detran::Material::SP_material   SP_material;
  typedef detran::Mesh::SP_mesh           SP_mesh;

  /*!
   *  \brief Constructor
   *  \param d
   */
  CartesianNodeDetran(const size_type  d,
                      const size_type  n,
                      const int        nodeid,
                      std::string      nodename,
                      const Point      nodeorigin,
                      vec2_size_type   so,
                      vec_size_type    po,
                      vec_size_type    ao,
                      vec_size_type    eo,
                      vec_dbl          nodewidth,
                      SP_db            nodedb,
                      SP_material      nodematerial,
                      SP_mesh          nodemesh);

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE
  //-------------------------------------------------------------------------//

  double color(Point point);

  //-------------------------------------------------------------------------//
  // DETRAN-SPECIFIC INTERFACE
  //-------------------------------------------------------------------------//

  SP_db db() const
  {
    return d_db;
  }

  SP_material material() const
  {
    return d_material;
  }

  SP_mesh mesh() const
  {
    return d_mesh;
  }

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Detran db
  SP_db d_db;

  /// Detran material
  SP_material d_material;

  /// Detran mesh
  SP_mesh d_mesh;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /// Default constructor is needed for serialization.
  CartesianNodeDetran(){}

#ifdef SERMENT_ENABLE_BOOST
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Base>(*this);
    ar & d_db;
    ar & d_material;
    ar & d_mesh;
  }
#endif

};

} // end namespace erme_geometry

#ifdef SERMENT_ENABLE_BOOST
BOOST_CLASS_EXPORT_KEY(erme_geometry::CartesianNodeDetran)
#endif

#endif // CARTESIANNODEDETRAN_HH_ 

//---------------------------------------------------------------------------//
//              end of file CartesianNodeDetran.hh
//---------------------------------------------------------------------------//
