// Copyright (C) 2011-2021 Vincent Heuveline
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

#include "alignedhexahedrontransformation.h"
#include "common/log.h"
#include <cmath>
#include <iomanip>

/// \author Philipp Gerstner

namespace hiflow {
namespace doffem {

// Reordering of vertices to make transformation coorespond to mesh
// ordering, with (0,0,0) mapped to vertex 0, and (1,1,1) mapped to vertex 7.

template < class DataType >
AlignedHexahedronTransformation< DataType >::AlignedHexahedronTransformation(
    int gdim)
    : CellTransformation< DataType >(gdim) {
  // define decomposition of hexahedron into 5 tetrahedrons
  //
  //        4 --------------- 7
  //       /|                /|
  //      / |               / |
  //     /  |z             /  |
  //    5 --------------- 6   |
  //    |   |             |   |
  //    |   |       y     |   |
  //    |   0 ------------|-- 3
  //    |  /              |  /
  //    | /x              | /
  //    |/                |/
  //    1 --------------- 2

  this->origin_ = 0;
  this->dx_ = 0.;
  this->dy_ = 0.;
  this->dz_ = 0.;
  this->x0_ = 0.;
  this->y0_ = 0.;
  this->z0_ = 0.;
#ifndef NDEBUG
  T3_ = new TriLinearHexahedronTransformation< DataType >(3);
#endif
}

template < class DataType >
std::vector< std::vector< DataType > >
AlignedHexahedronTransformation< DataType >::get_reference_coordinates() const {
  std::vector< std::vector< DataType > > ref(8);
  for (int i = 0; i < 8; ++i) {
    ref[i].resize(3, 0.);
  }
  ref[1][0] = 1.;

  ref[2][0] = 1.;
  ref[2][1] = 1.;

  ref[3][1] = 1.;

  ref[4][2] = 1.;

  ref[5][0] = 1.;
  ref[5][2] = 1.;

  ref[6][0] = 1.;
  ref[6][1] = 1.;
  ref[6][2] = 1.;

  ref[7][1] = 1.;
  ref[7][2] = 1.;

  return ref;
}

template < class DataType >
void AlignedHexahedronTransformation< DataType >::reinit(
    const Coord &coord_vtx) {
  this->coord_vtx_ = coord_vtx;
#ifdef SAVE_MODE
  std::vector< DataType > max_d(3, -1e12);
  std::vector< DataType > min_d(3, 1e12);
  std::vector< DataType > coord_sum(8, 0.);

  for (int p = 0; p < 8; ++p) {
    for (int d = 0; d < 3; ++d) {
      if (this->coord_vtx_[this->ij2ind(p, d)] > max_d[d]) {
        max_d[d] = this->coord_vtx_[this->ij2ind(p, d)];
      }
      if (this->coord_vtx_[this->ij2ind(p, d)] < min_d[d]) {
        min_d[d] = this->coord_vtx_[this->ij2ind(p, d)];
      }
      coord_sum[p] += this->coord_vtx_[this->ij2ind(p, d)];
    }
  }

  this->dx_ = max_d[0] - min_d[0];
  this->dy_ = max_d[1] - min_d[1];
  this->dz_ = max_d[2] - min_d[2];

  // find origin
  DataType min_sum = 1.e12;
  this->origin_ = 0;
  for (int p = 0; p < 8; ++p) {
    if (coord_sum[p] < min_sum) {
      min_sum = coord_sum[p];
      this->origin_ = p;
    }
  }
#else
  this->dx_ = this->coord_vtx_[this->ij2ind(3, 0)] -
              this->coord_vtx_[this->ij2ind(0, 0)];
  this->dy_ = this->coord_vtx_[this->ij2ind(1, 1)] -
              this->coord_vtx_[this->ij2ind(0, 1)];
  this->dz_ = this->coord_vtx_[this->ij2ind(4, 2)] -
              this->coord_vtx_[this->ij2ind(0, 2)];
  this->origin_ = 0;
#endif
  this->x0_ = this->coord_vtx_[this->ij2ind(this->origin_, 0)];
  this->y0_ = this->coord_vtx_[this->ij2ind(this->origin_, 1)];
  this->z0_ = this->coord_vtx_[this->ij2ind(this->origin_, 2)];

  assert(this->dx_ > 0.);
  assert(this->dy_ > 0.);
  assert(this->dz_ > 0.);
#ifndef NDEBUG
  this->T3_->reinit(coord_vtx);
#endif
}

template < class DataType >
bool AlignedHexahedronTransformation< DataType >::inverse(
    DataType x_phy, DataType y_phy, DataType &x_ref, DataType &y_ref) const {
  throw "This cell transformation does not support 2d inversion!\n";
  return false;
}

template < class DataType >
bool AlignedHexahedronTransformation< DataType >::inverse(
    DataType x_phy, DataType &x_ref) const {
  throw "This cell transformation does not support 1d inversion!\n";
  return false;
}

template < class DataType >
bool AlignedHexahedronTransformation< DataType >::inverse(
    DataType x_phy, DataType y_phy, DataType z_phy, DataType &x_ref,
    DataType &y_ref, DataType &z_ref) const {
  // Note (Philipp G): I know, y_ref = (x_phy ... and x_ref = (y_phy ... looks
  // wrong and the other way round should be more reasonable. However, there's a
  // problem somewhere in the DOF/FEM module related to checking which (local)
  // DOFs lie on a specific interface. Unfortunately, it seems that this
  // procedure does NOT take into account the underlying cell transformation.
  // Instead, this procedure assumes, for example, that those dofs that lie on
  // the reference facet with normal, say n = (0, 1, 0), do also lie on the
  // physical facet with normal n = (0, 1, 0), even if the underlying cell
  // transformation performs a rotation in the xy-plane. To make it even worse,
  // the DOF/FEM module assumes that a change in the x-axis in the reference
  // cell leads to a change in the y -axis in the physical cell. I guess this is
  // related to the internal dof ordering. To make a long story short: the
  // assumed DOF ordering is different from the assumed cell vertex ordering by
  // a rotation of 90 degress in the xy-plane.

  y_ref = (x_phy - this->x0_) / this->dx_;
  x_ref = (y_phy - this->y0_) / this->dy_;
  z_ref = (z_phy - this->z0_) / this->dz_;

#ifdef SAVE_MODE
  if (std::abs(x_ref) < this->eps_) {
    x_ref = 0.;
  } else if (std::abs(x_ref - 1.) < this->eps_) {
    x_ref = 1.;
  }

  if (std::abs(y_ref) < this->eps_) {
    y_ref = 0.;
  } else if (std::abs(y_ref - 1.) < this->eps_) {
    y_ref = 1.;
  }

  if (std::abs(z_ref) < this->eps_) {
    z_ref = 0.;
  } else if (std::abs(z_ref - 1.) < this->eps_) {
    z_ref = 1.;
  }
#endif
#ifndef NDEBUG
  Coord coord_ref(3);
  coord_ref[0] = x_ref;
  coord_ref[1] = y_ref;
  coord_ref[2] = z_ref;

  DataType x = this->x(coord_ref);
  DataType y = this->y(coord_ref);
  DataType z = this->z(coord_ref);

  assert(std::abs(x - x_phy) < 1e-12);
  assert(std::abs(y - y_phy) < 1e-12);
  assert(std::abs(z - z_phy) < 1e-12);
#if 0
            DataType x_ref3, y_ref3, z_ref3;
            bool found_std = this->T3_->inverse ( x_phy, y_phy, z_phy, x_ref3, y_ref3, z_ref3);
            DataType diff = std::abs ( x_ref - x_ref3 ) + std::abs ( y_ref - y_ref3 ) + std::abs ( z_ref - z_ref3 );
            if ( diff > 1e-13 && found_std )
            {
                LOG_DEBUG(0, "Aligned inverse : " << x_ref << ", " << y_ref << ", " << z_ref );
                LOG_DEBUG(0, "Standard inverse: " << x_ref3 << ", " << y_ref3 << ", " << z_ref3 );
            }
#endif
#endif
  return true;
}

template < class DataType >
DataType
AlignedHexahedronTransformation< DataType >::x(const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
#ifdef NDEBUG
  return this->x0_ + coord_ref[1] * this->dx_;
#else
  DataType x = this->x0_ + coord_ref[1] * this->dx_;
  DataType x3 = this->T3_->x(coord_ref);
  DataType diff = std::abs(x - x3);
  if (diff > 1e-13) {
    LOG_DEBUG(0, "Aligned x : " << x);
    LOG_DEBUG(0, "Standard x: " << x3);
  }
  return x;
#endif
}

template < class DataType >
DataType
AlignedHexahedronTransformation< DataType >::x_x(const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
#ifdef NDEBUG
  return 0.;
#else
  DataType x_x = 0.;
  DataType x_x3 = this->T3_->x_x(coord_ref);
  DataType diff = std::abs(x_x - x_x3);
  if (diff > 1e-13) {
    LOG_DEBUG(0, "Aligned x_x : " << x_x);
    LOG_DEBUG(0, "Standard x_x: " << x_x3);
  }
  return x_x;
#endif
}

template < class DataType >
DataType
AlignedHexahedronTransformation< DataType >::x_y(const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
#ifdef NDEBUG
  return this->dx_;
#else
  DataType x_y = this->dx_;
  DataType x_y3 = this->T3_->x_y(coord_ref);
  DataType diff = std::abs(x_y - x_y3);
  if (diff > 1e-13) {
    LOG_DEBUG(0, "Aligned x_y : " << x_y);
    LOG_DEBUG(0, "Standard x_y: " << x_y3);
  }
  return x_y;
#endif
}

template < class DataType >
DataType
AlignedHexahedronTransformation< DataType >::x_z(const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  return 0.;
}

template < class DataType >
DataType AlignedHexahedronTransformation< DataType >::x_xx(
    const Coord &coord_ref) const {
  return 0.;
}

template < class DataType >
DataType AlignedHexahedronTransformation< DataType >::x_xy(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  return 0.;
}

template < class DataType >
DataType AlignedHexahedronTransformation< DataType >::x_xz(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  return 0.;
}

template < class DataType >
DataType AlignedHexahedronTransformation< DataType >::x_yy(
    const Coord &coord_ref) const {
  return 0.;
}

template < class DataType >
DataType AlignedHexahedronTransformation< DataType >::x_yz(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  return 0.;
}

template < class DataType >
DataType AlignedHexahedronTransformation< DataType >::x_zz(
    const Coord &coord_ref) const {
  return 0.;
}

template < class DataType >
DataType
AlignedHexahedronTransformation< DataType >::y(const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  return this->y0_ + coord_ref[0] * this->dy_;
}

template < class DataType >
DataType
AlignedHexahedronTransformation< DataType >::y_x(const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  return this->dy_;
}

template < class DataType >
DataType
AlignedHexahedronTransformation< DataType >::y_y(const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  return 0.;
}

template < class DataType >
DataType
AlignedHexahedronTransformation< DataType >::y_z(const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  return 0.;
}

template < class DataType >
DataType AlignedHexahedronTransformation< DataType >::y_xx(
    const Coord &coord_ref) const {
  return 0.;
}

template < class DataType >
DataType AlignedHexahedronTransformation< DataType >::y_xy(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  return 0.;
}

template < class DataType >
DataType AlignedHexahedronTransformation< DataType >::y_xz(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  return 0.;
}

template < class DataType >
DataType AlignedHexahedronTransformation< DataType >::y_yy(
    const Coord &coord_ref) const {
  return 0.;
}

template < class DataType >
DataType AlignedHexahedronTransformation< DataType >::y_yz(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  return 0.;
}

template < class DataType >
DataType AlignedHexahedronTransformation< DataType >::y_zz(
    const Coord &coord_ref) const {
  return 0.;
}

template < class DataType >
DataType
AlignedHexahedronTransformation< DataType >::z(const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  return this->z0_ + coord_ref[2] * this->dz_;
}

template < class DataType >
DataType
AlignedHexahedronTransformation< DataType >::z_x(const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  return 0.;
}

template < class DataType >
DataType
AlignedHexahedronTransformation< DataType >::z_y(const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  return 0;
}

template < class DataType >
DataType
AlignedHexahedronTransformation< DataType >::z_z(const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  return this->dz_;
}

template < class DataType >
DataType AlignedHexahedronTransformation< DataType >::z_xx(
    const Coord &coord_ref) const {
  return 0.;
}

template < class DataType >
DataType AlignedHexahedronTransformation< DataType >::z_xy(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  return 0.;
}

template < class DataType >
DataType AlignedHexahedronTransformation< DataType >::z_xz(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  return 0.;
}

template < class DataType >
DataType AlignedHexahedronTransformation< DataType >::z_yy(
    const Coord &coord_ref) const {
  return 0.;
}

template < class DataType >
DataType AlignedHexahedronTransformation< DataType >::z_yz(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  return 0.;
}

template < class DataType >
DataType AlignedHexahedronTransformation< DataType >::z_zz(
    const Coord &coord_ref) const {
  return 0.;
}

template < class DataType >
bool AlignedHexahedronTransformation< DataType >::contains_reference_point(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == this->gdim_);
  return coord_ref[0] >= -this->eps_ && coord_ref[0] <= 1. + this->eps_ &&
         coord_ref[1] >= -this->eps_ && coord_ref[1] <= 1. + this->eps_ &&
         coord_ref[2] >= -this->eps_ && coord_ref[2] <= 1. + this->eps_;
}

template class AlignedHexahedronTransformation< double >;
template class AlignedHexahedronTransformation< float >;

} // namespace doffem
} // namespace hiflow
