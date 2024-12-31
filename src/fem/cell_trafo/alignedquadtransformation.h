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

#ifndef __FEM_ALIGNED_QUAD_TRANSFORMATION_H_
#define __FEM_ALIGNED_QUAD_TRANSFORMATION_H_

#include "fem/cell_trafo/bilinearquadtransformation.h"
#include "fem/cell_trafo/cell_transformation.h"
#include "mesh/geometric_tools.h"
#include "common/log.h"
#include "common/vector_algebra.h"
#include <cmath>
#include <iomanip>
#include <cassert>

namespace hiflow {
namespace doffem {

///
/// \class AlignedHexahedronTransformation alignedhexahedrontransformation.h
/// \brief Trilinear transformation mapping from reference to physical cell for
/// a Hexahedron which is assumed to be axis aligned \author Philipp Gerstner
///

template < class DataType, int DIM >
class AlignedQuadTransformation : public CellTransformation< DataType, DIM > {
public:
  typedef Vec<DIM, DataType> Coord;

  explicit AlignedQuadTransformation(ConstRefCellPtr<DataType, DIM> ref_cell);

  bool differs_by_translation_from (ConstCellTrafoPtr<DataType, DIM> rhs) const;
  
  ~ AlignedQuadTransformation()
  {
  }

#if 0
  void copy_from (CellTransformation<DataType, DIM> const * rhs)
  {
    CellTransformation<DataType, DIM>::copy_from(rhs);
    
    AlignedQuadTransformation<DataType, DIM> const * tmp_rhs = 
      dynamic_cast<AlignedQuadTransformation<DataType, DIM> const *>(rhs);
    
    assert (tmp_rhs != 0);
  }  
#endif

  void reinit(const std::vector<DataType> &coord_vtx);

  bool inverse(Coord co_phy, Coord &co_ref) const;

  DataType x(const Coord &coord_ref) const;
  DataType x_x(const Coord &coord_ref) const;
  DataType x_y(const Coord &coord_ref) const;

  DataType y(const Coord &coord_ref) const;
  DataType y_x(const Coord &coord_ref) const;
  DataType y_y(const Coord &coord_ref) const;

protected:
  Mat<DIM, DIM, DataType> A_;
  Mat<DIM, DIM, DataType> Ainv_; 
};

// Reordering of vertices to make transformation coorespond to mesh
// ordering, with (0,0,0) mapped to vertex 0, and (1,1,1) mapped to vertex 7.

template < class DataType, int DIM >
AlignedQuadTransformation< DataType, DIM >::AlignedQuadTransformation(ConstRefCellPtr<DataType, DIM> ref_cell)
    : CellTransformation< DataType, DIM >(ref_cell) 
{
  this->order_ = 1;
  this->fixed_ref_cell_type_ = REF_CELL_QUAD_STD;
  this->name_ = "AlignedQuad";
  
  //
  //        0 --------------- 3
  //       /        y        /
  //      /x                /
  //     /                 /
  //    1 --------------- 2

}

template < class DataType, int DIM >
bool AlignedQuadTransformation< DataType, DIM >::differs_by_translation_from (ConstCellTrafoPtr<DataType, DIM> rhs) const
{
  assert(this->ref_cell_);
  
  if (this->ref_cell_->type() != rhs->get_ref_cell()->type())
  {
    return false;
  }
  
  std::vector< Coord > ref_coords = this->ref_cell_->get_coords();
  assert (ref_coords.size() == 4);
  
  Coord my_p01 = this->transform(ref_coords[0]) - this->transform(ref_coords[1]); 
  Coord rhs_p01 = rhs->transform(ref_coords[0]) - rhs->transform(ref_coords[1]); 
  if (my_p01 != rhs_p01)
  {
    return false;
  }
  Coord my_p12 = this->transform(ref_coords[1]) - this->transform(ref_coords[2]); 
  Coord rhs_p12 = rhs->transform(ref_coords[1]) - rhs->transform(ref_coords[2]); 
  if (my_p12 != rhs_p12)
  {
    return false;
  }
  Coord my_p23 = this->transform(ref_coords[2]) - this->transform(ref_coords[3]); 
  Coord rhs_p23 = rhs->transform(ref_coords[2]) - rhs->transform(ref_coords[3]); 
  if (my_p23 != rhs_p23)
  {
    return false;
  }
  Coord my_p30 = this->transform(ref_coords[3]) - this->transform(ref_coords[0]); 
  Coord rhs_p30 = rhs->transform(ref_coords[3]) - rhs->transform(ref_coords[0]); 
  if (my_p30 != rhs_p30)
  {
    return false;
  }
  return true;
}

template < class DataType, int DIM >
void AlignedQuadTransformation< DataType, DIM >::reinit(const std::vector<DataType> &coord_vtx) {
  assert (DIM == 2);
  
  assert (coord_vtx.size() == DIM * 4);
  
  this->coord_vtx_.clear();
  this->coord_vtx_.resize(4);
  
  for (int i=0; i<4; ++i)
  {
    for (int d=0; d<2; ++d)
    {
      this->coord_vtx_[i][d] = coord_vtx[this->ij2ind(i,d)];
    }
  }

  //assert (mesh::is_aligned_rectangular_quad(coord_vtx));
  assert (mesh::is_parallelogram(coord_vtx));
  
  this->A_(0,0) = this->coord_vtx_[1][0] - this->coord_vtx_[0][0];
  this->A_(1,0) = this->coord_vtx_[1][1] - this->coord_vtx_[0][1];
  this->A_(0,1) = this->coord_vtx_[3][0] - this->coord_vtx_[0][0];
  this->A_(1,1) = this->coord_vtx_[3][1] - this->coord_vtx_[0][1];

#ifndef NDEBUG   
  DataType determ = this->A_(0,0) * this->A_(1,1) - this->A_(0,1) * this->A_(1,0);
  DataType sca = this->A_(0,0) * this->A_(0,1) + this->A_(1,0) * this->A_(1,1); 
  assert (std::abs(determ) > 1e-6);
  //assert (std::abs(sca) < 1e-10);
#endif

  inv(this->A_, this->Ainv_);
}

template < class DataType, int DIM >
bool AlignedQuadTransformation< DataType, DIM >::inverse(Coord co_phy, Coord &co_ref) const 
{
  assert (DIM == 2);

  Coord rhs = co_phy - this->coord_vtx_[0];
  this->Ainv_.VectorMult(rhs, co_ref);
  
  const DataType eps = this->ref_cell_->eps();
  if (std::abs(co_ref[0]) < eps) 
  {
    co_ref[0] = 0.;
  } 
  else if (std::abs(co_ref[0] - 1.) < eps) 
  {
    co_ref[0] = 1.;
  }

  if (std::abs(co_ref[1]) < eps) 
  {
    co_ref[1] = 0.;
  } 
  else if (std::abs(co_ref[1] - 1.) < eps) 
  {
    co_ref[1] = 1.;
  }

#ifndef NDEBUG
  DataType x = this->x(co_ref);
  DataType y = this->y(co_ref);

  assert(std::abs(x - co_phy[0]) < 1e-12);
  assert(std::abs(y - co_phy[1]) < 1e-12);
#endif
  return true;
}

template < class DataType, int DIM >
DataType AlignedQuadTransformation< DataType, DIM >::x(const Coord &coord_ref) const {
  assert (DIM == 2);
  return this->coord_vtx_[0][0] + this->A_(0,0) * coord_ref[0] + this->A_(0,1) * coord_ref[1]; 
}

template < class DataType, int DIM >
DataType AlignedQuadTransformation< DataType, DIM >::x_x(const Coord &coord_ref) const {
  assert (DIM == 2);
  return this->A_(0,0);
}

template < class DataType, int DIM >
DataType AlignedQuadTransformation< DataType, DIM >::x_y(const Coord &coord_ref) const {
  assert (DIM == 2);
  return this->A_(0,1);
}

template < class DataType, int DIM >
DataType AlignedQuadTransformation< DataType, DIM >::y(const Coord &coord_ref) const {
  assert (DIM == 2);
  return this->coord_vtx_[0][1] + this->A_(1,0) * coord_ref[0] + this->A_(1,1) * coord_ref[1]; 
}

template < class DataType, int DIM >
DataType AlignedQuadTransformation< DataType, DIM >::y_x(const Coord &coord_ref) const {
  assert (DIM == 2);
  return this->A_(1,0); 
}

template < class DataType, int DIM >
DataType AlignedQuadTransformation< DataType, DIM >::y_y(const Coord &coord_ref) const {
  assert (DIM == 2);
  return this->A_(1,1); 
}

} // namespace doffem
} // namespace hiflow

#endif
