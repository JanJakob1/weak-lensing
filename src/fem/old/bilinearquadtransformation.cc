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

#include "bilinearquadtransformation.h"
#include "common/log.h"
#include "lineartriangletransformation.h"
#include <cmath>
#include <iomanip>

/// \author Michael Schick<br>Martin Baumann<br>Simon Gawlok<br>Philipp Gerstner

namespace hiflow {
namespace doffem {

template < class DataType >
BiLinearQuadTransformation< DataType >::BiLinearQuadTransformation(int gdim)
    : CellTransformation< DataType >(gdim) {
  // define decomposition of quadrilateral into 2 triangle
  //
  //        0 --------------- 3
  //       /        y        /
  //      /x                /
  //     /                 /
  //    1 --------------- 2

  // triangle 0:
  //        0 --------------- 3
  //       /        y
  //      /x
  //     /
  //    1

  triangle_decomp_ind_[0][0] = 0;
  triangle_decomp_ind_[0][1] = 1;
  triangle_decomp_ind_[0][2] = 3;

  // triangle 1:
  //                          3
  //                         /
  //                       x/
  //            y          /
  //    1 --------------- 2

  triangle_decomp_ind_[1][0] = 2;
  triangle_decomp_ind_[1][1] = 3;
  triangle_decomp_ind_[1][2] = 1;
}

template < class DataType >
std::vector< std::vector< DataType > >
BiLinearQuadTransformation< DataType >::get_reference_coordinates() const {
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
bool BiLinearQuadTransformation< DataType >::inverse_by_decomposition(
    DataType x_phy, DataType y_phy, DataType &x_ref, DataType &y_ref) const {
  LinearTriangleTransformation< DataType > tria_trafo(2);
  DataType x_ref_tri, y_ref_tri;

  bool found_by_decomposition = false;
  int tria_ind = -1;

  // loop through tetrahedron decomposition
  for (int t = 0; t < 2; ++t) {
    // build linear tetrahedron transformation
    Coord tria_coord(6, 0.);
    for (int v = 0; v < 3; ++v) {
      for (int d = 0; d < 2; ++d) {
        tria_coord[this->ij2ind(v, d)] =
            this->coord_vtx_[this->ij2ind(this->triangle_decomp_ind_[t][v], d)];
      }
    }
    tria_trafo.reinit(tria_coord);

    // compute reference coordinates w.r.t. tetrahedron
    bool tria_success = tria_trafo.inverse(x_phy, y_phy, x_ref_tri, y_ref_tri);
    if (!tria_success) {
      continue;
    }

    // check whether reference coordinates are contained in tetrahedron
    Coord ref_coord(2);
    ref_coord[0] = x_ref_tri;
    ref_coord[1] = y_ref_tri;

    // if yes: convert reference coordinates
    if (tria_trafo.contains_reference_point(ref_coord)) {
      switch (t) {
      case 0:
        x_ref = x_ref_tri;
        y_ref = y_ref_tri;
        break;
      case 1:
        x_ref = 1. - x_ref_tri;
        y_ref = 1. - y_ref_tri;
        break;
      }
      found_by_decomposition = true;
      tria_ind = t;
      break;
    }
  }
  // check if reference point found by decomposition is sufficiently accurate
  Coord ref_coord(2);
  ref_coord[0] = x_ref;
  ref_coord[1] = y_ref;

  DataType x_phy_test = this->x(ref_coord);
  DataType y_phy_test = this->y(ref_coord);

  DataType residual = std::sqrt((x_phy_test - x_phy) * (x_phy_test - x_phy) +
                                (y_phy_test - y_phy) * (y_phy_test - y_phy));

#ifndef NDEBUG
  LOG_DEBUG(2, "==============");
  LOG_DEBUG(3, "Physical point " << x_phy << " " << y_phy
                                 << " is contained in triangle " << tria_ind);

  if (tria_ind == -1) {
    LOG_DEBUG(2, "Physical point " << x_phy << " " << y_phy
                                   << " is contained in triangle " << tria_ind);
    for (int v = 0; v < 4; ++v) {
      LOG_DEBUG(2, "vertex " << v << " : "
                             << this->coord_vtx_[this->ij2ind(v, 0)] << ", "
                             << this->coord_vtx_[this->ij2ind(v, 1)]);
    }
  }

  LOG_DEBUG(
      3, "Inverse cell-trafo based on triangle decomposition yields residual = "
             << residual);
#endif

  if (residual < 5. * std::numeric_limits< DataType >::epsilon()) {
    return true;
  } else {
    return false;
  }
}

template < class DataType >
bool BiLinearQuadTransformation< DataType >::inverse(DataType x_phy,
                                                     DataType y_phy,
                                                     DataType &x_ref,
                                                     DataType &y_ref) const {
  bool success = this->inverse_by_decomposition(x_phy, y_phy, x_ref, y_ref);

  if (success) {
    // reference point obtained by decomposition is sufficiently accurate ->
    // exit
    return true;
  } else {
    // reference point obtained by decomposition is not accurate enough -> use
    // as initial value for Newton method
    return this->inverse_newton_2d(x_phy, y_phy, x_ref, y_ref, x_ref, y_ref);
    /*
    // Newton method with old initial value
    if (!success)
    {
        success = this->inverse_newton_2d ( x_phy, y_phy, x_ref, y_ref );
    }
    * */
  }
}

template < class DataType >
bool BiLinearQuadTransformation< DataType >::inverse(
    DataType x_phy, DataType y_phy, DataType z_phy, DataType &x_ref,
    DataType &y_ref, DataType &z_ref) const {
  throw "This cell transformation does not support 3d inversion!\n";
  return false;
}

template < class DataType >
bool BiLinearQuadTransformation< DataType >::inverse(DataType x_phy,
                                                     DataType &x_ref) const {
  throw "This cell transformation does not support 1d inversion!\n";
  return false;
}

template < class DataType >
DataType
BiLinearQuadTransformation< DataType >::x(const Coord &coord_ref) const {
  assert(coord_ref.size() >= 2);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_1 = coord_ref[1];

  return this->coord_vtx_[this->ij2ind(0, 0)] * (1. - coord_0) *
             (1. - coord_1) +
         this->coord_vtx_[this->ij2ind(1, 0)] * coord_0 * (1. - coord_1) +
         this->coord_vtx_[this->ij2ind(2, 0)] * coord_0 * coord_1 +
         this->coord_vtx_[this->ij2ind(3, 0)] * (1. - coord_0) * coord_1;
}

template < class DataType >
DataType
BiLinearQuadTransformation< DataType >::x_x(const Coord &coord_ref) const {
  assert(coord_ref.size() >= 2);
  const DataType coord_1 = coord_ref[1];

  return this->coord_vtx_[this->ij2ind(0, 0)] * (-1. + coord_1) +
         this->coord_vtx_[this->ij2ind(1, 0)] * (1. - coord_1) +
         this->coord_vtx_[this->ij2ind(2, 0)] * coord_1 -
         this->coord_vtx_[this->ij2ind(3, 0)] * coord_1;
}

template < class DataType >
DataType
BiLinearQuadTransformation< DataType >::x_y(const Coord &coord_ref) const {
  assert(coord_ref.size() >= 2);
  const DataType coord_0 = coord_ref[0];

  return this->coord_vtx_[this->ij2ind(0, 0)] * (-1. + coord_0) -
         this->coord_vtx_[this->ij2ind(1, 0)] * coord_0 +
         this->coord_vtx_[this->ij2ind(2, 0)] * coord_0 +
         this->coord_vtx_[this->ij2ind(3, 0)] * (1. - coord_0);
}

template < class DataType >
DataType
BiLinearQuadTransformation< DataType >::x_xx(const Coord &coord_ref) const {
  return 0.;
}

template < class DataType >
DataType
BiLinearQuadTransformation< DataType >::x_xy(const Coord &coord_ref) const {
  assert(coord_ref.size() >= 2);
  return this->coord_vtx_[this->ij2ind(0, 0)] -
         this->coord_vtx_[this->ij2ind(1, 0)] +
         this->coord_vtx_[this->ij2ind(2, 0)] -
         this->coord_vtx_[this->ij2ind(3, 0)];
}

template < class DataType >
DataType
BiLinearQuadTransformation< DataType >::x_yy(const Coord &coord_ref) const {
  return 0.;
}

template < class DataType >
DataType
BiLinearQuadTransformation< DataType >::y(const Coord &coord_ref) const {
  assert(coord_ref.size() >= 2);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_1 = coord_ref[1];

  return this->coord_vtx_[this->ij2ind(0, 1)] * (1. - coord_0) *
             (1. - coord_1) +
         this->coord_vtx_[this->ij2ind(1, 1)] * coord_0 * (1. - coord_1) +
         this->coord_vtx_[this->ij2ind(2, 1)] * coord_0 * coord_1 +
         this->coord_vtx_[this->ij2ind(3, 1)] * (1. - coord_0) * coord_1;
}

template < class DataType >
DataType
BiLinearQuadTransformation< DataType >::y_x(const Coord &coord_ref) const {
  assert(coord_ref.size() >= 2);
  const DataType coord_1 = coord_ref[1];

  return this->coord_vtx_[this->ij2ind(0, 1)] * (-1. + coord_1) +
         this->coord_vtx_[this->ij2ind(1, 1)] * (1. - coord_1) +
         this->coord_vtx_[this->ij2ind(2, 1)] * coord_1 -
         this->coord_vtx_[this->ij2ind(3, 1)] * coord_1;
}

template < class DataType >
DataType
BiLinearQuadTransformation< DataType >::y_y(const Coord &coord_ref) const {
  assert(coord_ref.size() >= 2);
  const DataType coord_0 = coord_ref[0];

  return this->coord_vtx_[this->ij2ind(0, 1)] * (-1. + coord_0) -
         this->coord_vtx_[this->ij2ind(1, 1)] * coord_0 +
         this->coord_vtx_[this->ij2ind(2, 1)] * coord_0 +
         this->coord_vtx_[this->ij2ind(3, 1)] * (1. - coord_0);
}

template < class DataType >
DataType
BiLinearQuadTransformation< DataType >::y_xx(const Coord &coord_ref) const {
  return 0.;
}

template < class DataType >
DataType
BiLinearQuadTransformation< DataType >::y_xy(const Coord &coord_ref) const {
  assert(coord_ref.size() >= 2);

  return this->coord_vtx_[this->ij2ind(0, 1)] -
         this->coord_vtx_[this->ij2ind(1, 1)] +
         this->coord_vtx_[this->ij2ind(2, 1)] -
         this->coord_vtx_[this->ij2ind(3, 1)];
}

template < class DataType >
DataType
BiLinearQuadTransformation< DataType >::y_yy(const Coord &coord_ref) const {
  return 0.;
}

// TODO avoid trunc

template < class DataType >
bool BiLinearQuadTransformation< DataType >::contains_reference_point(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == this->gdim_);
  return coord_ref[0] >= -this->eps_ && coord_ref[0] <= 1. + this->eps_ &&
         coord_ref[1] >= -this->eps_ && coord_ref[1] <= 1. + this->eps_;
}

template class BiLinearQuadTransformation< double >;
template class BiLinearQuadTransformation< float >;

} // namespace doffem
} // namespace hiflow
