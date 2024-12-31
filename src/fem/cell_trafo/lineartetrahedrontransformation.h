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

#ifndef __FEM_LINEAR_TETRAHEDRON_TRANSFORMATION_H_
#define __FEM_LINEAR_TETRAHEDRON_TRANSFORMATION_H_

#include "fem/cell_trafo/cell_transformation.h"
#include <cmath>
#include <iomanip>

namespace hiflow {
namespace doffem {

///
/// \class LinearTetrahedronTransformation lineartetrahedrontransformation.h
/// \brief Linear transformation mapping from reference to physical cell for a
/// Tetrahedron \author Michael Schick<br>Martin Baumann
///

template < class DataType, int DIM >
class LinearTetrahedronTransformation : public CellTransformation< DataType, DIM > {
public:
  typedef Vec<DIM, DataType> Coord;

  explicit LinearTetrahedronTransformation(ConstRefCellPtr<DataType, DIM> ref_cell);

  bool differs_by_translation_from (ConstCellTrafoPtr<DataType, DIM> rhs) const;
    
  bool inverse(Coord co_phy, Coord &co_ref) const;

  void reinit(const std::vector<DataType> &coord_vtx);
  
  DataType x(const Coord &coord_ref) const;
  DataType x_x(const Coord &coord_ref) const;
  DataType x_y(const Coord &coord_ref) const;
  DataType x_z(const Coord &coord_ref) const;

  DataType y(const Coord &coord_ref) const;
  DataType y_x(const Coord &coord_ref) const;
  DataType y_y(const Coord &coord_ref) const;
  DataType y_z(const Coord &coord_ref) const;

  DataType z(const Coord &coord_ref) const;
  DataType z_x(const Coord &coord_ref) const;
  DataType z_y(const Coord &coord_ref) const;
  DataType z_z(const Coord &coord_ref) const;

private:
  Mat<DIM, DIM, DataType> A_;
  Mat<DIM, DIM, DataType> Ainv_;
};

template < class DataType, int DIM >
LinearTetrahedronTransformation< DataType, DIM >::LinearTetrahedronTransformation(ConstRefCellPtr<DataType, DIM> ref_cell)
    : CellTransformation< DataType, DIM >(ref_cell) 
{
  this->order_ = 1;
  this->fixed_ref_cell_type_ = REF_CELL_TET_STD; 
  this->name_ = "Tet";
}

template < class DataType, int DIM >
bool LinearTetrahedronTransformation< DataType, DIM >::differs_by_translation_from (ConstCellTrafoPtr<DataType, DIM> rhs) const
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
  Coord my_p02 = this->transform(ref_coords[0]) - this->transform(ref_coords[2]); 
  Coord rhs_p02 = rhs->transform(ref_coords[0]) - rhs->transform(ref_coords[2]); 
  if (my_p02 != rhs_p02)
  {
    return false;
  }
  Coord my_p03 = this->transform(ref_coords[0]) - this->transform(ref_coords[3]); 
  Coord rhs_p03 = rhs->transform(ref_coords[0]) - rhs->transform(ref_coords[3]); 
  if (my_p03 != rhs_p03)
  {
    return false;
  }
  return true;
}

template < class DataType, int DIM >
void LinearTetrahedronTransformation< DataType, DIM >::reinit(const std::vector<DataType> &coord_vtx) {
  assert (DIM == 3);
  
  assert (coord_vtx.size() == DIM * 4);
  
  this->coord_vtx_.clear();
  this->coord_vtx_.resize(4);
  
  for (int i=0; i<4; ++i)
  {
    for (int d=0; d<3; ++d)
    {
      this->coord_vtx_[i][d] = coord_vtx[this->ij2ind(i,d)];
    }
  }
  
  this->A_(0,0) = this->coord_vtx_[1][0] - this->coord_vtx_[0][0];
  this->A_(1,0) = this->coord_vtx_[1][1] - this->coord_vtx_[0][1];
  this->A_(2,0) = this->coord_vtx_[1][2] - this->coord_vtx_[0][2];
  this->A_(0,1) = this->coord_vtx_[2][0] - this->coord_vtx_[0][0];
  this->A_(1,1) = this->coord_vtx_[2][1] - this->coord_vtx_[0][1];
  this->A_(2,1) = this->coord_vtx_[2][2] - this->coord_vtx_[0][2];
  this->A_(0,2) = this->coord_vtx_[3][0] - this->coord_vtx_[0][0];
  this->A_(1,2) = this->coord_vtx_[3][1] - this->coord_vtx_[0][1];
  this->A_(2,2) = this->coord_vtx_[3][2] - this->coord_vtx_[0][2];

#ifndef NDEBUG
  DataType determ = det(this->A_);
  assert (std::abs(determ) > 1e-6);
#endif

  inv(this->A_, this->Ainv_);
}

template < class DataType, int DIM >
bool LinearTetrahedronTransformation< DataType, DIM >::inverse(Coord co_phy, Coord &co_ref) const {
  assert (DIM == 3);
  
  Coord rhs = co_phy - this->coord_vtx_[0];
  this->Ainv_.VectorMult(rhs, co_ref);

  const DataType eps = this->ref_cell_->eps();
  for (int d=0; d<DIM; ++d)
  {
    if (std::abs(co_ref[d]) < eps) 
    {
      co_ref[d] = 0.;
    } 
    else if (std::abs(co_ref[d] - 1.) < eps) 
    {
      co_ref[d] = 1.;
    }
  }
  
/*
  DataType a11 = this->coord_vtx_[1][0] - this->coord_vtx_[0][0];
  DataType a12 = this->coord_vtx_[2][0] - this->coord_vtx_[0][0];
  DataType a13 = this->coord_vtx_[3][0] - this->coord_vtx_[0][0];

  DataType a21 = this->coord_vtx_[1][1] - this->coord_vtx_[0][1];
  DataType a22 = this->coord_vtx_[2][1] - this->coord_vtx_[0][1];
  DataType a23 = this->coord_vtx_[3][1] - this->coord_vtx_[0][1];

  DataType a31 = this->coord_vtx_[1][2] - this->coord_vtx_[0][2];
  DataType a32 = this->coord_vtx_[2][2] - this->coord_vtx_[0][2];
  DataType a33 = this->coord_vtx_[3][2] - this->coord_vtx_[0][2];

  DataType det = a11 * (a33 * a22 - a32 * a23) - a21 * (a33 * a12 - a32 * a13) + a31 * (a23 * a12 - a22 * a13);

  assert(det != 0.0);

  co_ref[0] = (1.0 / det) * ((a33 * a22 - a32 * a23) *
                             (co_phy[0] - this->coord_vtx_[0][0]) -
                         (a33 * a12 - a32 * a13) *
                             (co_phy[1] - this->coord_vtx_[0][1]) +
                         (a23 * a12 - a22 * a13) *
                             (co_phy[2] - this->coord_vtx_[0][2]));

  co_ref[1] = (1.0 / det) * (-(a33 * a21 - a31 * a23) *
                             (co_phy[0] - this->coord_vtx_[0][0]) +
                         (a33 * a11 - a31 * a13) *
                             (co_phy[1] - this->coord_vtx_[0][1]) -
                         (a23 * a11 - a21 * a13) *
                             (co_phy[2] - this->coord_vtx_[0][2]));

  co_ref[2] = (1.0 / det) * ((a32 * a21 - a31 * a22) *
                             (co_phy[0] - this->coord_vtx_[0][0]) -
                         (a32 * a11 - a31 * a12) *
                             (co_phy[1] - this->coord_vtx_[0][1]) +
                         (a22 * a11 - a21 * a12) *
                             (co_phy[2] - this->coord_vtx_[0][2]));
*/
  return this->contains_reference_point(co_ref);
}

template < class DataType, int DIM >
DataType LinearTetrahedronTransformation< DataType, DIM >::x(const Coord &coord_ref) const {
  assert (DIM == 3);
  return (this->coord_vtx_[1][0] - this->coord_vtx_[0][0]) * coord_ref[0] +
         (this->coord_vtx_[2][0] - this->coord_vtx_[0][0]) * coord_ref[1] +
         (this->coord_vtx_[3][0] - this->coord_vtx_[0][0]) * coord_ref[2] +
          this->coord_vtx_[0][0];
}

template < class DataType, int DIM >
DataType LinearTetrahedronTransformation< DataType, DIM >::x_x(const Coord &coord_ref) const {
  assert (DIM == 3);
  return (this->coord_vtx_[1][0] - this->coord_vtx_[0][0]);
}

template < class DataType, int DIM >
DataType LinearTetrahedronTransformation< DataType, DIM >::x_y(const Coord &coord_ref) const {
  assert (DIM == 3);
  return (this->coord_vtx_[2][0] - this->coord_vtx_[0][0]);
}

template < class DataType, int DIM >
DataType LinearTetrahedronTransformation< DataType, DIM >::x_z(const Coord &coord_ref) const {
  assert (DIM == 3);
  return (this->coord_vtx_[3][0] - this->coord_vtx_[0][0]);
}

template < class DataType, int DIM >
DataType LinearTetrahedronTransformation< DataType, DIM >::y(const Coord &coord_ref) const {
  assert (DIM == 3);
  return (this->coord_vtx_[1][1] - this->coord_vtx_[0][1]) * coord_ref[0] +
         (this->coord_vtx_[2][1] - this->coord_vtx_[0][1]) * coord_ref[1] +
         (this->coord_vtx_[3][1] - this->coord_vtx_[0][1]) * coord_ref[2] +
         this->coord_vtx_[0][1];
}

template < class DataType, int DIM >
DataType LinearTetrahedronTransformation< DataType, DIM >::y_x(const Coord &coord_ref) const {
  assert (DIM == 3);
  return (this->coord_vtx_[1][1] - this->coord_vtx_[0][1]);
}

template < class DataType, int DIM >
DataType LinearTetrahedronTransformation< DataType, DIM >::y_y(const Coord &coord_ref) const {
  assert (DIM == 3);
  return (this->coord_vtx_[2][1] - this->coord_vtx_[0][1]);
}

template < class DataType, int DIM >
DataType LinearTetrahedronTransformation< DataType, DIM >::y_z(const Coord &coord_ref) const {
  assert (DIM == 3);
  return (this->coord_vtx_[3][1] - this->coord_vtx_[0][1]);
}

template < class DataType, int DIM >
DataType LinearTetrahedronTransformation< DataType, DIM >::z(const Coord &coord_ref) const {
  assert (DIM == 3);
  return (this->coord_vtx_[1][2] - this->coord_vtx_[0][2]) * coord_ref[0] +
         (this->coord_vtx_[2][2] - this->coord_vtx_[0][2]) * coord_ref[1] +
         (this->coord_vtx_[3][2] - this->coord_vtx_[0][2]) * coord_ref[2] +
         this->coord_vtx_[0][2];
}

template < class DataType, int DIM >
DataType LinearTetrahedronTransformation< DataType, DIM >::z_x(const Coord &coord_ref) const {
  assert (DIM == 3);
  return (this->coord_vtx_[1][2] - this->coord_vtx_[0][2]);
}

template < class DataType, int DIM >
DataType LinearTetrahedronTransformation< DataType, DIM >::z_y(const Coord &coord_ref) const {
  assert (DIM == 3);
  return (this->coord_vtx_[2][2] - this->coord_vtx_[0][2]);
}

template < class DataType, int DIM >
DataType LinearTetrahedronTransformation< DataType, DIM >::z_z(const Coord &coord_ref) const {
  assert (DIM == 3);
  return (this->coord_vtx_[3][2] - this->coord_vtx_[0][2]);
}

} // namespace doffem
} // namespace hiflow

#endif
