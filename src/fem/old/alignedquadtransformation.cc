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

#include "alignedquadtransformation.h"
#include "common/log.h"
#include <cmath>
#include <iomanip>

/// \author Philipp Gerstner

namespace hiflow {
namespace doffem {

// Reordering of vertices to make transformation coorespond to mesh
// ordering, with (0,0,0) mapped to vertex 0, and (1,1,1) mapped to vertex 7.

template < class DataType >
AlignedQuadTransformation< DataType >::AlignedQuadTransformation(int gdim)
    : CellTransformation< DataType >(gdim) {
  //
  //        0 --------------- 3
  //       /        y        /
  //      /x                /
  //     /                 /
  //    1 --------------- 2

  this->origin_ = 0;
  this->dx_ = 0.;
  this->dy_ = 0.;
  this->x0_ = 0.;
  this->y0_ = 0.;

#ifndef NDEBUG
  T2_ = new BiLinearQuadTransformation< DataType >(2);
#endif
}

template < class DataType >
std::vector< std::vector< DataType > >
AlignedQuadTransformation< DataType >::get_reference_coordinates() const {
  std::vector< std::vector< DataType > > ref(4);
  for (int i = 0; i < 4; ++i) {
    ref[i].resize(2, 0.);
  }
  ref[1][0] = 1.;

  ref[2][0] = 1.;
  ref[2][1] = 1.;

  ref[3][1] = 1.;

  return ref;
}

template < class DataType >
void AlignedQuadTransformation< DataType >::reinit(const Coord &coord_vtx) {
  this->coord_vtx_ = coord_vtx;

#ifdef SAVE_MODE
  std::vector< DataType > max_d(2, -1e12);
  std::vector< DataType > min_d(2, 1e12);
  std::vector< DataType > coord_sum(4, 0.);

  for (int p = 0; p < 4; ++p) {
    for (int d = 0; d < 2; ++d) {
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

  // find origin
  DataType min_sum = 1.e12;
  this->origin_ = 0;
  for (int p = 0; p < 4; ++p) {
    if (coord_sum[p] < min_sum) {
      min_sum = coord_sum[p];
      this->origin_ = p;
    }
  }
#else
  this->dx_ = this->coord_vtx_[this->ij2ind(1, 0)] -
              this->coord_vtx_[this->ij2ind(0, 0)];
  this->dy_ = this->coord_vtx_[this->ij2ind(3, 1)] -
              this->coord_vtx_[this->ij2ind(0, 1)];
  this->origin_ = 0;
#endif
  this->x0_ = this->coord_vtx_[this->ij2ind(this->origin_, 0)];
  this->y0_ = this->coord_vtx_[this->ij2ind(this->origin_, 1)];

  assert(this->dx_ > 0.);
  assert(this->dy_ > 0.);
#ifndef NDEBUG
  this->T2_->reinit(coord_vtx);
#endif
}

template < class DataType >
bool AlignedQuadTransformation< DataType >::inverse(DataType x_phy,
                                                    DataType y_phy,
                                                    DataType &x_ref,
                                                    DataType &y_ref) const {
  x_ref = (x_phy - this->x0_) / this->dx_;
  y_ref = (y_phy - this->y0_) / this->dy_;

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
#endif
#ifndef NDEBUG
  Coord coord_ref(3);
  coord_ref[0] = x_ref;
  coord_ref[1] = y_ref;

  DataType x = this->x(coord_ref);
  DataType y = this->y(coord_ref);

  assert(std::abs(x - x_phy) < 1e-12);
  assert(std::abs(y - y_phy) < 1e-12);
#if 0
            DataType x_ref2, y_ref2;
            bool found_std = this->T2_->inverse ( x_phy, y_phy, x_ref2, y_ref2);
            DataType diff = std::abs ( x_ref - x_ref2 ) + std::abs ( y_ref - y_ref2 );
            if ( diff > 1e-14 && found_std )
            {
                LOG_DEBUG(1, "Aligned inverse : " << x_ref << ", " << y_ref  );
                LOG_DEBUG(1, "Standard inverse: " << x_ref2 << ", " << y_ref2 );
            }
#endif
#endif
  return true;
}

template < class DataType >
bool AlignedQuadTransformation< DataType >::inverse(DataType x_phy,
                                                    DataType &x_ref) const {
  throw "This cell transformation does not support 1d inversion!\n";
  return false;
}

template < class DataType >
bool AlignedQuadTransformation< DataType >::inverse(
    DataType x_phy, DataType y_phy, DataType z_phy, DataType &x_ref,
    DataType &y_ref, DataType &z_ref) const {
  throw "This cell transformation does not support 3d inversion!\n";
  return false;
}

template < class DataType >
DataType
AlignedQuadTransformation< DataType >::x(const Coord &coord_ref) const {
  assert(coord_ref.size() >= 2);
#ifdef NDEBUG
  return this->x0_ + coord_ref[0] * this->dx_;
#else
  DataType x = this->x0_ + coord_ref[0] * this->dx_;
  DataType x2 = this->T2_->x(coord_ref);
  DataType diff = std::abs(x - x2);
  if (diff > 1e-14) {
    LOG_DEBUG(1, "Aligned x : " << x);
    LOG_DEBUG(1, "Standard x: " << x2);
  }
  return x;
#endif
}

template < class DataType >
DataType
AlignedQuadTransformation< DataType >::x_x(const Coord &coord_ref) const {
  assert(coord_ref.size() >= 2);
#ifdef NDEBUG
  return this->dx_;
#else
  DataType x_x = this->dx_;
  DataType x_x2 = this->T2_->x_x(coord_ref);
  DataType diff = std::abs(x_x - x_x2);
  if (diff > 1e-14) {
    LOG_DEBUG(1, "Aligned x_x : " << x_x);
    LOG_DEBUG(1, "Standard x_x: " << x_x2);
  }
  return x_x;
#endif
}

template < class DataType >
DataType
AlignedQuadTransformation< DataType >::x_y(const Coord &coord_ref) const {
  assert(coord_ref.size() >= 2);
#ifdef NDEBUG
  return 0.;
#else
  DataType x_y = 0.;
  DataType x_y2 = this->T2_->x_y(coord_ref);
  DataType diff = std::abs(x_y - x_y2);
  if (diff > 1e-14) {
    LOG_DEBUG(1, "Aligned x_y : " << x_y);
    LOG_DEBUG(1, "Standard x_y: " << x_y2);
  }
  return x_y;
#endif
}

template < class DataType >
DataType
AlignedQuadTransformation< DataType >::x_xx(const Coord &coord_ref) const {
  return 0.;
}

template < class DataType >
DataType
AlignedQuadTransformation< DataType >::x_xy(const Coord &coord_ref) const {
  assert(coord_ref.size() >= 2);
  return 0.;
}

template < class DataType >
DataType
AlignedQuadTransformation< DataType >::x_yy(const Coord &coord_ref) const {
  return 0.;
}

template < class DataType >
DataType
AlignedQuadTransformation< DataType >::y(const Coord &coord_ref) const {
  assert(coord_ref.size() >= 2);
  return this->y0_ + coord_ref[1] * this->dy_;
}

template < class DataType >
DataType
AlignedQuadTransformation< DataType >::y_x(const Coord &coord_ref) const {
  assert(coord_ref.size() >= 2);
  return 0.;
}

template < class DataType >
DataType
AlignedQuadTransformation< DataType >::y_y(const Coord &coord_ref) const {
  assert(coord_ref.size() >= 2);
  return this->dy_;
}

template < class DataType >
DataType
AlignedQuadTransformation< DataType >::y_xx(const Coord &coord_ref) const {
  return 0.;
}

template < class DataType >
DataType
AlignedQuadTransformation< DataType >::y_xy(const Coord &coord_ref) const {
  assert(coord_ref.size() >= 2);
  return 0.;
}

template < class DataType >
DataType
AlignedQuadTransformation< DataType >::y_yy(const Coord &coord_ref) const {
  return 0.;
}

template < class DataType >
bool AlignedQuadTransformation< DataType >::contains_reference_point(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == this->gdim_);
  return coord_ref[0] >= -this->eps_ && coord_ref[0] <= 1. + this->eps_ &&
         coord_ref[1] >= -this->eps_ && coord_ref[1] <= 1. + this->eps_;
}

template class AlignedQuadTransformation< double >;
template class AlignedQuadTransformation< float >;

} // namespace doffem
} // namespace hiflow
