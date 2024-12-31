// Copyright (C) 2011-2017 Vincent Heuveline
//
// HiFlow3 is free software: you can redistribute it and/or modify it under the
// terms of the European Union Public Licence (EUPL) v1.2 as published by the
// European Union or (at your option) any later version.
//
// HiFlow3 is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the European Union Public Licence (EUPL) v1.2 for
// more details.
//
// You should have received a copy of the European Union Public Licence (EUPL)
// v1.2 along with HiFlow3.  If not, see
// <https://joinup.ec.europa.eu/page/eupl-text-11-12>.

#ifndef __FEM_REF_CELL_H_
#define __FEM_REF_CELL_H_

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include <cmath>

#include "common/macros.h"
#include "common/vector_algebra.h"
#include "common/log.h"
#include "mesh/cell_type.h"
#include "dof/dof_fem_types.h"

namespace hiflow {
namespace doffem {

///
/// \class ReferenceCell reference_cell.h
/// \brief Ancestor class of different reference cells used for defining a FE
/// \author Philipp Gerstner


template < class DataType, int DIM > 
class RefCell 
{
public:
  typedef Vec<DIM, DataType> Coord;

  /// Default Constructor
  RefCell()
  : topo_cell_(nullptr), 
    type_(REF_CELL_NOT_SET),
    // std::numeric_limits<double>::epsilon ( ) = 2.22045e-16
    // 1e3 is used, because same tolerance is used in function compute_weights(),
    // see numbering_lagrange.cc
    eps_ (1.e3 * std::numeric_limits< DataType >::epsilon()) 
  {}

  /// Default Destructor  
  // TODO: delete cell object?
  virtual ~RefCell() 
  {}

  /// Tpological dimension of cell
  inline size_t tdim() const
  {
    assert (this->topo_cell_ != nullptr);
    return this->topo_cell_->tdim();
  }
  
  inline RefCellType type () const
  {
    return this->type_;
  }

  inline mesh::CellType const * topo_cell () const
  {
    return this->topo_cell_;
  }

  inline mesh::CellType::Tag tag () const
  {
    assert (this->topo_cell_ != nullptr);
    return this->topo_cell_->tag();
  }
  
  inline size_t num_vertices () const 
  {
    return this->coord_vertices_.size();
  }

  inline DataType eps() const
  {
    return this->eps_;
  }

  inline size_t num_entities(int tdim) const 
  {
    return this->topo_cell_->num_entities(tdim);
  }
  
  virtual std::vector< Coord > get_coords () const = 0;

  /// check whether coordinates of underlying reference cell coincide with 
  /// given set of points
  bool ref_coord_match( const std::vector<Coord> & test_coord) const
  {
    std::vector<Coord> my_coord = this->get_coords();
    if (my_coord.size() == 0 || test_coord.size() == 0)
    {
      return false;
    }
    if (my_coord.size() != test_coord.size() )
    {
      return false;
    }
    for (size_t p = 0; p<test_coord.size(); ++p)
    {
      if (norm (my_coord[p] - test_coord[p]) > this->eps_)
      {
        return false;
      }
    }
    return true;
  }

  virtual bool contains_point ( const Coord & pt) const = 0;

  virtual void compute_facet_projection_matrix (int facet_number, Mat< DIM, DIM - 1, DataType >& proj) const = 0;

  virtual void compute_facet_normal (int facet_number, Vec< DIM, DataType >& n) const = 0;

  virtual void compute_facet_tangents (int facet_number, Vec< DIM, DataType >& t1, Vec< DIM, DataType >& t2) const = 0;

  /// Map point on reference facet (line, triangle, quad) to coordinates on facet on reference cell with given facet_nr
  virtual Coord Tf (size_t facet_nr, const Vec<DIM-1,DataType> & ref_pt) const = 0;

  virtual DataType ds (size_t facet_nr, const Vec<DIM-1,DataType> & ref_pt) const = 0;

  virtual std::string get_quad_name_cell_gauss (bool use_economical) const = 0;
  virtual std::string get_quad_name_facet_gauss (int facet_nr, bool use_economical) const = 0;
  virtual mesh::CellType::Tag facet_tag(int facet_nr) const = 0;

protected:
  /// Storing an instance of the reference cell
  mesh::CellType const * topo_cell_;

  RefCellType type_;

  /// tolerance for checking whether reference coordinates lie in cell
  DataType eps_;

};

template < class DataType, int DIM > 
class RefCellLineStd : public RefCell<DataType, DIM>
{
public:
  typedef Vec<DIM, DataType> Coord;
  
  /// Default Constructor
  RefCellLineStd()
  {
    this->topo_cell_ = &(mesh::CellType::get_instance(mesh::CellType::LINE));
    this->type_ = REF_CELL_LINE_STD;
  }

  std::vector< Coord > get_coords () const
  {
    std::vector<Coord> ref_coord(2);
    ref_coord[1][0] = 1.;

    return ref_coord;  
  }

  void compute_facet_projection_matrix (int facet_number, Mat< DIM, DIM - 1, DataType >& proj) const
  {
    assert (0);
  }

  void compute_facet_normal (int facet_number, Vec< DIM, DataType >& n) const
  {
    assert (0);
  }

  void compute_facet_tangents (int facet_number, Vec< DIM, DataType >& t1, Vec< DIM, DataType >& t2) const
  {
    assert (0);
  }

  bool contains_point (const Coord &pt) const 
  {
    return pt[0] >= -this->eps_ && pt[0] <= 1. + this->eps_;
  }

  Coord Tf (size_t facet_nr, const Vec<DIM-1,DataType> & ref_pt) const
  {
    assert (0);
    Coord tmp;
    return tmp;
  }

  DataType ds (size_t facet_nr, const Vec<DIM-1,DataType> & ref_pt) const
  {
    assert (0);
    return 0.;
  }

  std::string get_quad_name_cell_gauss (bool use_economical) const 
  {
    return "GaussLine";  
  }

  std::string get_quad_name_facet_gauss (int facet_nr, bool use_economical) const 
  {
    return "";  
  }

  mesh::CellType::Tag facet_tag(int facet_nr) const
  {
    return mesh::CellType::POINT;
  }

  /// Default Destructor
  virtual ~RefCellLineStd()
  {
  }
};

template < class DataType, int DIM > 
class RefCellQuadStd : public RefCell<DataType, DIM>
{
public:
  typedef Vec<DIM, DataType> Coord;
  
  /// Default Constructor
  RefCellQuadStd()
  {
    this->topo_cell_ = &(mesh::CellType::get_instance(mesh::CellType::QUADRILATERAL));
    this->type_ = REF_CELL_QUAD_STD;
  }

  std::vector< Coord > get_coords () const
  {
    std::vector<Coord> ref_coord(4);
    ref_coord[1][0] = 1.;
    ref_coord[2][0] = 1.;
    ref_coord[2][1] = 1.;
    ref_coord[3][1] = 1.;

    return ref_coord;  
  }

  void compute_facet_projection_matrix (int facet_number, Mat< DIM, DIM - 1, DataType >& proj) const
  {
    assert (DIM > 1);
    proj.Zeros();
    switch (facet_number) 
    {
      case 0: // bottom edge
      case 2: // top edge
        proj(0, 0) = 1.;
        break;
      case 1: // right edge
      case 3: // left edge
        proj(1, 0) = 1.;
        break;
      default:
        assert(0);
    }
  }

  void compute_facet_normal (int facet_number, Vec< DIM, DataType >& n) const
  {
    assert (DIM > 1);
    n.Zeros();
    switch (facet_number) 
    {
      case 0: // bottom edge
        n[1] = -1.;
        break;
      case 1: // right edge
        n[0] = 1.;
        break;
      case 2: // top edge
        n[1] = 1.;
        break;
      case 3: // left edge
        n[0] = -1.;
        break;
      default:
        assert(0);
    }
  }

  void compute_facet_tangents (int facet_number, Vec< DIM, DataType >& t1, Vec< DIM, DataType >& t2) const
  {
    t1.Zeros();
    t2.Zeros();
    NOT_YET_IMPLEMENTED;
  }

  Coord Tf (size_t facet_nr, const Vec<DIM-1,DataType> & ref_pt) const
  {
    assert (DIM == 2);
    Coord mapped_pt;
    switch (facet_nr)
    {
      case 0:
        mapped_pt[0] = ref_pt[0];
        break;
      case 1:
        mapped_pt[0] = 1.;
        mapped_pt[1] = ref_pt[0];
        break;
      case 2:
        mapped_pt[0] = ref_pt[0];
        mapped_pt[1] = 1.;
        break;
      case 3:
        mapped_pt[1] = ref_pt[0];
        break;
      default:
        assert(0);  
        break;
    }
    return mapped_pt;
  }

  DataType ds (size_t facet_nr, const Vec<DIM-1,DataType> & ref_pt) const
  {
    return 1.;
  }

  bool contains_point (const Coord &pt) const 
  {
    return pt[0] >= -this->eps_ && pt[0] <= 1. + this->eps_ &&
           pt[1] >= -this->eps_ && pt[1] <= 1. + this->eps_;
  }

  std::string get_quad_name_cell_gauss (bool use_economical) const 
  {
    if (use_economical)
    {
      return "EconomicalGaussQuadrilateral";
    }
    else
    {
      return "GaussQuadrilateral";
    }  
  }

  std::string get_quad_name_facet_gauss (int facet_nr, bool use_economical) const 
  {
    return "GaussLine";  
  }

  mesh::CellType::Tag facet_tag(int facet_nr) const
  {
    return mesh::CellType::LINE;
  }

  /// Default Destructor
  virtual ~RefCellQuadStd()
  {
  }
};

template < class DataType, int DIM > 
class RefCellTriStd : public RefCell<DataType, DIM>
{
public:
  typedef Vec<DIM, DataType> Coord;
  
  /// Default Constructor
  RefCellTriStd()
  {
    this->topo_cell_ = &(mesh::CellType::get_instance(mesh::CellType::TRIANGLE));
    this->type_ = REF_CELL_TRI_STD;
  }

  std::vector< Coord > get_coords () const
  {
    std::vector<Coord> ref_coord(3);
    ref_coord[1][0] = 1.;
    ref_coord[2][1] = 1.;

    return ref_coord;  
  }

  void compute_facet_projection_matrix (int facet_number, Mat< DIM, DIM - 1, DataType >& proj) const
  {
    assert (DIM > 1);
    proj.Zeros();
    switch (facet_number) 
    {
      case 0: // bottom edge
        proj(0, 0) = 1.;
        break;
      case 1: // diagonal
        proj(0, 0) = 1.;
        proj(1, 0) = -1.;
        break;
      case 2: // left edge
        proj(1, 0) = 1.;
        break;
      default:
        assert(0);
    }
  }

  void compute_facet_normal (int facet_number, Vec< DIM, DataType >& n) const
  {
    assert (DIM > 1);
    n.Zeros();
    switch (facet_number) 
    {
      case 0: // bottom edge
        n[1] = -1.;
        break;
      case 1: // diagonal
        n[0] = 1. / std::sqrt(2.);
        n[1] = 1. / std::sqrt(2.);
        break;
      case 2: // left edge
        n[0] = -1.;
        break;
      default:
        assert(0);
    }
  }

  void compute_facet_tangents (int facet_number, Vec< DIM, DataType >& t1, Vec< DIM, DataType >& t2) const
  {
    t1.Zeros();
    t2.Zeros();
    NOT_YET_IMPLEMENTED;
  }

  Coord Tf (size_t facet_nr, const Vec<DIM-1,DataType> & ref_pt) const
  {
    assert (DIM == 2);
    Coord mapped_pt;
    switch (facet_nr)
    {
      case 0:
        mapped_pt[0] = ref_pt[0];
        break;
      case 1:
        mapped_pt[0] = 1. - ref_pt[0];
        mapped_pt[1] = ref_pt[0];
        break;
      case 2:
        mapped_pt[1] = 1. - ref_pt[0];
        break;
      default:
        assert(0);  
        break;
    }
    return mapped_pt;
  }

  DataType ds (size_t facet_nr, const Vec<DIM-1,DataType> & ref_pt) const
  {
    switch (facet_nr)
    {
      case 1:
        return std::sqrt(2.);
        break;
      default:
        return 1.;
        break;
    }
  }

  bool contains_point (const Coord &pt) const 
  {
    return pt[0] >= -this->eps_ && pt[0] <= 1. + this->eps_ &&
           pt[1] >= -this->eps_ &&
           pt[1] <= 1. - pt[0] + this->eps_;
  }

  std::string get_quad_name_cell_gauss (bool use_economical) const 
  {
    return "GaussTriangle";
  }

  std::string get_quad_name_facet_gauss (int facet_nr, bool use_economical) const 
  {
    return "GaussLine";  
  }

  mesh::CellType::Tag facet_tag(int facet_nr) const
  {
    return mesh::CellType::LINE;
  }

  /// Default Destructor
  virtual ~RefCellTriStd()
  {
  }
};

template < class DataType, int DIM > 
class RefCellHexStd : public RefCell<DataType, DIM>
{
  
  //              f5 
  //        4 --------------- 7
  //       /|                /|
  //      / |       f4      / |
  //     /  |z,k           /  |
  //    5 --------------- 6   |
  //    |   |             |   |  f3
  // f2 |   |       y,j   |   |
  //    |   0 ------------|-- 3
  //    |  /              |  /
  //    | /x,i     f1     | /
  //    |/                |/
  //    1 --------------- 2
  //             f0

public:
  typedef Vec<DIM, DataType> Coord;
  
  /// Default Constructor
  RefCellHexStd()
  {
    this->topo_cell_ = &(mesh::CellType::get_instance(mesh::CellType::HEXAHEDRON));
    this->type_ = REF_CELL_HEX_STD;

  }

  /// Default Destructor
  virtual ~RefCellHexStd()
  {
  }

  std::vector< Coord > get_coords () const
  {
    std::vector<Coord> ref_coord(8);

    ref_coord[1][0] = 1.;
    ref_coord[2][0] = 1.;
    ref_coord[2][1] = 1.;
    ref_coord[3][1] = 1.;
    ref_coord[4][2] = 1.;
    ref_coord[5][0] = 1.;
    ref_coord[5][2] = 1.;
    ref_coord[6][0] = 1.;
    ref_coord[6][1] = 1.;
    ref_coord[6][2] = 1.;
    ref_coord[7][1] = 1.;
    ref_coord[7][2] = 1.;

    return ref_coord;  
  }

  void compute_facet_projection_matrix (int facet_number, Mat< DIM, DIM - 1, DataType >& proj) const
  {
    assert (DIM == 3);
    proj.Zeros();
    switch (facet_number) 
    {
      case 0: // bottom
      case 5: // top
        proj(0, 0) = 1.;
        proj(1, 0) = 0.;
        proj(2, 0) = 0.;

        proj(0, 1) = 0.;
        proj(1, 1) = 1.;
        proj(2, 1) = 0.;
        break;

      case 1: // front 1
      case 4: // back  4
        proj(0, 0) = 1.;
        proj(1, 0) = 0.;
        proj(2, 0) = 0.;

        proj(0, 1) = 0.;
        proj(1, 1) = 0.;
        proj(2, 1) = 1.;
        break;

      case 2: // left  2
      case 3: // right 3
        proj(0, 0) = 0.;
        proj(1, 0) = 1.;
        proj(2, 0) = 0.;

        proj(0, 1) = 0.;
        proj(1, 1) = 0.;
        proj(2, 1) = 1.;
        break;

      default:
        assert(0);
    }
  }

  void compute_facet_normal (int facet_number, Vec< DIM, DataType >& n) const
  {
    assert (DIM == 3);
    n.Zeros();
    switch (facet_number) 
    {
    case 0: // bottom
      n[2] = -1.;
      break;
    case 5: // top
      n[2] = 1.;
      break;
    case 1: // front
      n[1] = -1.;
      break;
    case 4: // back
      n[1] = 1.;
      break;
    case 2: // left
      n[0] = -1.;
      break;
    case 3: // right
      n[0] = 1.;
      break;
    default:
      assert(0);
    }
  }

  void compute_facet_tangents (int facet_number, Vec< DIM, DataType >& t1, Vec< DIM, DataType >& t2) const
  {
    t1.Zeros();
    t2.Zeros();
    NOT_YET_IMPLEMENTED;
  }

  Coord Tf (size_t facet_nr, const Vec<DIM-1,DataType> & ref_pt) const
  {
    assert (DIM == 3);
    Coord mapped_pt;
    switch (facet_nr)
    {
      case 0:
        mapped_pt[0] = ref_pt[0];
        mapped_pt[1] = ref_pt[1];
        break;
      case 1:
        mapped_pt[0] = 1.;
        mapped_pt[1] = ref_pt[0];
        mapped_pt[2] = ref_pt[1];
        break;
      case 2:
        mapped_pt[0] = ref_pt[0];
        mapped_pt[2] = ref_pt[1];
        break;
      case 3:
        mapped_pt[0] = ref_pt[0];
        mapped_pt[1] = 1.;
        mapped_pt[2] = ref_pt[1];
        break;
      case 4:
        mapped_pt[1] = ref_pt[0];
        mapped_pt[2] = ref_pt[1];
        break;
      case 5:
        mapped_pt[0] = ref_pt[0];
        mapped_pt[1] = ref_pt[1];
        mapped_pt[2] = 1.;
        break;
      default:
        assert(0);  
        break;
    }
    return mapped_pt;
  }

  DataType ds (size_t facet_nr, const Vec<DIM-1,DataType> & ref_pt) const
  {
    return 1.;
  }

  bool contains_point (const Coord &pt) const 
  {
    return pt[0] >= -this->eps_ && pt[0] <= 1. + this->eps_ &&
           pt[1] >= -this->eps_ && pt[1] <= 1. + this->eps_ &&
           pt[2] >= -this->eps_ && pt[2] <= 1. + this->eps_;
  }

  std::string get_quad_name_cell_gauss (bool use_economical) const 
  {
    if (use_economical)
    {
      return "EconomicalGaussHexahedron";
    }
    else
    {
      return "GaussHexahedron";
    }  
  }

  std::string get_quad_name_facet_gauss (int facet_nr, bool use_economical) const 
  {
    if (use_economical)
    {
      return "EconomicalGaussQuadrilateral";
    }
    else
    {
      return "GaussQuadrilateral";
    }    
  }

  mesh::CellType::Tag facet_tag(int facet_nr) const
  {
    return mesh::CellType::QUADRILATERAL;
  }

};

template < class DataType, int DIM > 
class RefCellTetStd : public RefCell<DataType, DIM>
{
    typedef Vec<DIM, DataType> Coord;
    
public:
//     z^
//      |\
//      | \
//      |  \
//      |   \
//   f1 | f2 \
//      |_____\ ->y
//     /
//    /  f0
//   /
//  x

  /// Default Constructor
  RefCellTetStd()
  {
    this->topo_cell_ = &(mesh::CellType::get_instance(mesh::CellType::TETRAHEDRON));
    this->type_ = REF_CELL_TET_STD;
  }

  /// Default Destructor
  virtual ~RefCellTetStd()
  {
  }

  std::vector< Coord > get_coords () const
  {
    std::vector<Coord> ref_coord(4);

    ref_coord[1][0] = 1.;
    ref_coord[2][1] = 1.;
    ref_coord[3][2] = 1.;

    return ref_coord;  
  }

  void compute_facet_projection_matrix (int facet_number, Mat< DIM, DIM - 1, DataType >& proj) const
  {
    assert (DIM == 3);
    proj.Zeros();
    switch (facet_number) 
    {
      case 0: // bottom
        proj(0, 0) = 1.;
        proj(1, 0) = 0.;
        proj(2, 0) = 0.;

        proj(0, 1) = 0.;
        proj(1, 1) = 1.;
        proj(2, 1) = 0.;
        break;

      case 1: // front
        proj(0, 0) = 1.;
        proj(1, 0) = 0.;
        proj(2, 0) = 0.;

        proj(0, 1) = 0.;
        proj(1, 1) = 0.;
        proj(2, 1) = 1.;
        break;

      case 2: // left
        proj(0, 0) = 0.;
        proj(1, 0) = 1.;
        proj(2, 0) = 0.;

        proj(0, 1) = 0.;
        proj(1, 1) = 0.;
        proj(2, 1) = 1.;
        break;

      case 3: // back
        /*proj(0, 0) = 1.;
        proj(1, 0) = -1.;
        proj(2, 0) = 0.;

        proj(0, 1) = 0.;
        proj(1, 1) = 1.;
        proj(2, 1) = -1.;*/
        proj(0, 0) = 1.;
        proj(1, 0) = 0.;
        proj(2, 0) = -1.;

        proj(0, 1) = 0.;
        proj(1, 1) = 1.;
        proj(2, 1) = -1.;
        break;

      default:
        assert(0);
    }
  }

  void compute_facet_normal (int facet_number, Vec< DIM, DataType >& n) const
  {
    assert (DIM == 3);
    n.Zeros();
    switch (facet_number) 
    {
      case 0: // bottom
        n[2] = -1.;
        break;
      case 1: // front
        n[1] = -1.;
        break;
      case 2: // left
        n[0] = -1.;
        break;
      case 3: // back
        n[0] = 1. / std::sqrt(3.);
        n[1] = 1. / std::sqrt(3.);
        n[2] = 1. / std::sqrt(3.);
        break;
      default:
        assert(0);
    }
  }

  void compute_facet_tangents (int facet_number, Vec< DIM, DataType >& t1, Vec< DIM, DataType >& t2) const
  {
    t1.Zeros();
    t2.Zeros();
    NOT_YET_IMPLEMENTED;
  }

  Coord Tf (size_t facet_nr, const Vec<DIM-1,DataType> & ref_pt) const
  {
    assert (DIM == 3);
    Coord mapped_pt;
    switch (facet_nr)
    {
      case 0:
        mapped_pt[0] = ref_pt[0];
        mapped_pt[1] = ref_pt[1];
        break;
      case 1:
        mapped_pt[0] = ref_pt[0];
        mapped_pt[2] = ref_pt[1];
        break;
      case 2:
        mapped_pt[1] = ref_pt[0];
        mapped_pt[2] = ref_pt[1];
        break;
      case 3:
        mapped_pt[0] = ref_pt[0];
        mapped_pt[1] = ref_pt[1];
        mapped_pt[2] = 1. - ref_pt[0] - ref_pt[1];
        break;
      default:
        assert(0);  
        break;
    }
    return mapped_pt;
  }

  DataType ds (size_t facet_nr, const Vec<DIM-1,DataType> & ref_pt) const
  {
    switch (facet_nr)
    {
      case 3:
        return std::sqrt(3.);
        break;
      default:
        return 1.;
        break;
    }
  }

  bool contains_point (const Coord &pt) const 
  {
    return pt[0] >= -this->eps_ && pt[0] <= 1. + this->eps_ &&
           pt[1] >= -this->eps_ &&
           pt[1] <= 1. - pt[0] + this->eps_ &&
           pt[2] >= -this->eps_ &&
           pt[2] <= 1. - pt[0] - pt[1] + this->eps_;
  }

  std::string get_quad_name_cell_gauss (bool use_economical) const 
  {
    return "GaussTetrahedron";
  }

  std::string get_quad_name_facet_gauss (int facet_nr, bool use_economical) const 
  {
    return "GaussTriangle";  
  }

  mesh::CellType::Tag facet_tag(int facet_nr) const
  {
    return mesh::CellType::TRIANGLE;
  }
};

template < class DataType, int DIM > 
class RefCellPyrStd : public RefCell<DataType, DIM>
{
public:
  typedef Vec<DIM, DataType> Coord;
  
  /// Default Constructor
  RefCellPyrStd()
  {
    this->topo_cell_ = &(mesh::CellType::get_instance(mesh::CellType::PYRAMID));
    this->type_ = REF_CELL_PYR_STD;
  }

  /// Default Destructor
  virtual ~RefCellPyrStd()
  {
  }

  std::vector< Coord > get_coords () const
  {
    std::vector<Coord> ref_coord(5);

    ref_coord[1][0] = 1.;
    ref_coord[2][0] = 1.;
    ref_coord[2][1] = 1.;
    ref_coord[3][1] = 1.;
    ref_coord[4][0] = 0.5;
    ref_coord[4][1] = 0.5;
    ref_coord[4][2] = 1.;

    return ref_coord;  
  }

  void compute_facet_projection_matrix (int facet_number, Mat< DIM, DIM - 1, DataType >& proj) const
  {
    assert (DIM ==3);
    proj.Zeros();
    switch (facet_number) 
    {
      case 0: // bottom
        proj(0, 0) = 1.;
        proj(1, 0) = 0.;
        proj(2, 0) = 0.;

        proj(1, 1) = 0.;
        proj(1, 1) = 1.;
        proj(2, 1) = 0.;
        break;

      case 1: // front
        proj(0, 0) = 1.;
        proj(1, 0) = 0.;
        proj(2, 0) = 0.;

        proj(0, 1) = 0.5;
        proj(1, 1) = 0.5;
        proj(2, 1) = 1.;
        break;

      case 2: // right
        proj(0, 0) = 0.;
        proj(1, 0) = 1.;
        proj(2, 0) = 0.;

        proj(0, 1) = -0.5;
        proj(1, 1) = 0.5;
        proj(2, 1) = 1.;
        break;

      case 3: // back
        proj(0, 0) = -1.;
        proj(1, 0) = 0.;
        proj(2, 0) = 0.;

        proj(0, 1) = -0.5;
        proj(1, 1) = -0.;
        proj(2, 1) = 1.;
        break;

      case 4: // left
        proj(0, 0) = 0.;
        proj(1, 0) = -1.;
        proj(2, 0) = 0.;

        proj(0, 1) = 0.5;
        proj(1, 1) = -0.5;
        proj(2, 1) = 1.;
        break;

      default:
        assert(0);
    }
  }

  void compute_facet_normal (int facet_number, Vec< DIM, DataType >& n) const
  {
    assert (DIM == 3);
    n.Zeros();

    switch (facet_number) 
    {
      case 0: // bottom
        n[2] = -1.;
        break;
      case 1: // front
        n[1] = -2. / std::sqrt(3.);
        n[2] = 1. / std::sqrt(3.);
        break;
      case 2: // right
        n[0] = 2. / std::sqrt(3.);
        n[2] = 1. / std::sqrt(3.);
        break;
      case 3: // back
        n[1] = 2. / std::sqrt(3.);
        n[2] = 1. / std::sqrt(3.);
        break;
      case 4: // left
        n[0] = -2. / std::sqrt(3.);
        n[2] = 1. / std::sqrt(3.);
        break;
      default:
        assert(0);
    }
  }

  void compute_facet_tangents (int facet_number, Vec< DIM, DataType >& t1, Vec< DIM, DataType >& t2) const
  {
    t1.Zeros();
    t2.Zeros();

    NOT_YET_IMPLEMENTED;
  }

  Coord Tf (size_t facet_nr, const Vec<DIM-1,DataType> & ref_pt) const
  {
    NOT_YET_IMPLEMENTED;
    Coord mapped_pt;
    return mapped_pt;
  }

  DataType ds (size_t facet_nr, const Vec<DIM-1,DataType> & ref_pt) const
  {
    NOT_YET_IMPLEMENTED;
    return 0.;
  }

  bool contains_point (const Coord &pt) const 
  {
    return pt[0] >= + 0.5 * pt[2] - this->eps_
       &&  pt[0] <= 1. - 0.5 * pt[2] + this->eps_
       &&  pt[1] >= 0. + 0.5 * pt[2] - this->eps_
       &&  pt[1] <= 1. - 0.5 * pt[2] + this->eps_
       &&  pt[2] >= -this->eps_ && pt[2] <= 1. + this->eps_;
  }
  
  std::string get_quad_name_cell_gauss (bool use_economical) const 
  {
    return "GaussPyramid";
  }

  std::string get_quad_name_facet_gauss (int facet_nr, bool use_economical) const 
  {
    if (facet_nr == 0) 
    {
      if (use_economical)
      {
        return "EconomicalGaussQuadrilateral";
      }
      else
      {
        return "GaussQuadrilateral";
      }
    }
    else 
    {
      return "GaussTriangle";
    }  
  }

  mesh::CellType::Tag facet_tag(int facet_nr) const
  {
    if (facet_nr == 0)  
    {
      return mesh::CellType::QUADRILATERAL;
    }
    else
    {
      return mesh::CellType::TRIANGLE;
    }
  }
};


} // namespace doffem
} // namespace hiflow
#endif
