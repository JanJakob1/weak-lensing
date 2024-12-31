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

#include "trilinearhexahedrontransformation.h"
#include "common/log.h"
#include "lineartetrahedrontransformation.h"
#include <cmath>
#include <iomanip>

/// \author Michael Schick<br>Martin Baumann<br>Simon Gawlok<br>Philipp Gerstner

namespace hiflow {
namespace doffem {

// Reordering of vertices to make transformation coorespond to mesh
// ordering, with (0,0,0) mapped to vertex 0, and (1,1,1) mapped to vertex 7.

template < class DataType >
TriLinearHexahedronTransformation<
    DataType >::TriLinearHexahedronTransformation(int gdim)
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

  // tetrahedron 0:
  //                      6
  //                      |
  //                      |
  //                     z|   3
  //                      |  /
  //                      | /x
  //             y        |/
  //    1 --------------- 2

  tetrahedron_decomp_ind_[0][0] = 2;
  tetrahedron_decomp_ind_[0][1] = 3;
  tetrahedron_decomp_ind_[0][2] = 1;
  tetrahedron_decomp_ind_[0][3] = 6;

  // tetrahedron 1:
  //        4
  //        |
  //        |
  //        |z
  //        |
  //        |
  //        |
  //        0 --------------- 3
  //       /        y
  //     x/
  //     /
  //    1

  tetrahedron_decomp_ind_[1][0] = 0;
  tetrahedron_decomp_ind_[1][1] = 1;
  tetrahedron_decomp_ind_[1][2] = 3;
  tetrahedron_decomp_ind_[1][3] = 4;

  // tetrahedron 2:
  //        4
  //       /
  //     x/
  //     /       y
  //    5 --------------- 6
  //    |
  //    |
  //    |z
  //    |
  //    |
  //    |
  //    1

  tetrahedron_decomp_ind_[2][0] = 5;
  tetrahedron_decomp_ind_[2][1] = 4;
  tetrahedron_decomp_ind_[2][2] = 6;
  tetrahedron_decomp_ind_[2][3] = 1;

  // tetrahedron 3:
  //        4 --------------- 7
  //                y        /|
  //                       x/ |
  //                       /  |
  //                      6   |z
  //                          |
  //                          |
  //                          3

  tetrahedron_decomp_ind_[3][0] = 7;
  tetrahedron_decomp_ind_[3][1] = 6;
  tetrahedron_decomp_ind_[3][2] = 4;
  tetrahedron_decomp_ind_[3][3] = 3;

  // tetrahedron 4:
  //        4 -
  //             -
  //               z-
  //                   -
  //                      6
  //                   -    \x
  //                -        \
            //              -           3
  //            -y
  //         -
  //       -
  //    1

  tetrahedron_decomp_ind_[4][0] = 6;
  tetrahedron_decomp_ind_[4][1] = 3;
  tetrahedron_decomp_ind_[4][2] = 1;
  tetrahedron_decomp_ind_[4][3] = 4;
}

template < class DataType >
std::vector< std::vector< DataType > >
TriLinearHexahedronTransformation< DataType >::get_reference_coordinates()
    const {
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
bool TriLinearHexahedronTransformation< DataType >::inverse_by_decomposition(
    DataType x_phy, DataType y_phy, DataType z_phy, DataType &x_ref,
    DataType &y_ref, DataType &z_ref) const {
  LinearTetrahedronTransformation< DataType > tetra_trafo(3);
  DataType x_ref_tet, y_ref_tet, z_ref_tet;

  bool found_by_decomposition = false;
  int tetra_ind = -1;
  // loop through tetrahedron decomposition
  for (int t = 0; t < 5; ++t) {
    // build linear tetrahedron transformation
    Coord tetra_coord(12, 0.);
    for (int v = 0; v < 4; ++v) {
      for (int d = 0; d < 3; ++d) {
        tetra_coord[this->ij2ind(v, d)] = this->coord_vtx_[this->ij2ind(
            this->tetrahedron_decomp_ind_[t][v], d)];
      }
    }
    tetra_trafo.reinit(tetra_coord);

    // compute reference coordinates w.r.t. tetrahedron
    bool tetra_success = tetra_trafo.inverse(x_phy, y_phy, z_phy, x_ref_tet,
                                             y_ref_tet, z_ref_tet);
    if (!tetra_success) {
      continue;
    }

    // check whether reference coordinates are contained in tetrahedron
    Coord ref_coord(3);
    ref_coord[0] = x_ref_tet;
    ref_coord[1] = y_ref_tet;
    ref_coord[2] = z_ref_tet;

    // if yes: convert reference coordinates
    if (tetra_trafo.contains_reference_point(ref_coord)) {
      switch (t) {
      case 0:
        x_ref = 1. - x_ref_tet;
        y_ref = 1. - y_ref_tet;
        z_ref = z_ref_tet;
        break;
      case 1:
        x_ref = x_ref_tet;
        y_ref = y_ref_tet;
        z_ref = z_ref_tet;
        break;
      case 2:
        x_ref = 1. - x_ref_tet;
        y_ref = y_ref_tet;
        z_ref = 1. - z_ref_tet;
        break;
      case 3:
        x_ref = x_ref_tet;
        y_ref = 1. - y_ref_tet;
        z_ref = 1. - z_ref_tet;
        break;
      case 4:
        x_ref = 1. - x_ref_tet - z_ref_tet;
        y_ref = 1. - y_ref_tet - z_ref_tet;
        z_ref = 1. - x_ref_tet - y_ref_tet;
        break;
      }
      found_by_decomposition = true;
      tetra_ind = t;
      break;
    }
  }
  // check if reference point found by decomposition is sufficiently accurate
  Coord ref_coord(3);
  ref_coord[0] = x_ref;
  ref_coord[1] = y_ref;
  ref_coord[2] = z_ref;

  DataType x_phy_test = this->x(ref_coord);
  DataType y_phy_test = this->y(ref_coord);
  DataType z_phy_test = this->z(ref_coord);

  DataType residual = std::sqrt((x_phy_test - x_phy) * (x_phy_test - x_phy) +
                                (y_phy_test - y_phy) * (y_phy_test - y_phy) +
                                (z_phy_test - z_phy) * (z_phy_test - z_phy));

#ifndef NDEBUG
  LOG_DEBUG(3, "==============");
  LOG_DEBUG(3, "Physical point " << x_phy << " " << y_phy << " " << z_phy
                                 << " is contained in tetrahedron "
                                 << tetra_ind);

  if (tetra_ind == -1) {
    LOG_DEBUG(3, "Physical point " << x_phy << " " << y_phy << " " << z_phy
                                   << " is not contained in any tetrahedron ");
    for (int v = 0; v < 8; ++v) {
      LOG_DEBUG(3, "vertex " << v << " : "
                             << this->coord_vtx_[this->ij2ind(v, 0)] << ", "
                             << this->coord_vtx_[this->ij2ind(v, 1)] << ", "
                             << this->coord_vtx_[this->ij2ind(v, 2)]);
    }
  }

  LOG_DEBUG(
      3,
      "Inverse cell-trafo based on tetrahedron decomposition yields residual = "
          << residual);
#endif

  if (residual < 5. * std::numeric_limits< DataType >::epsilon()) {
    return true;
  } else {
    return false;
  }
}

template < class DataType >
bool TriLinearHexahedronTransformation< DataType >::inverse(
    DataType x_phy, DataType y_phy, DataType &x_ref, DataType &y_ref) const {
  throw "This cell transformation does not support 2d inversion!\n";
  return false;
}

template < class DataType >
bool TriLinearHexahedronTransformation< DataType >::inverse(
    DataType x_phy, DataType &x_ref) const {
  throw "This cell transformation does not support 1d inversion!\n";
  return false;
}

template < class DataType >
bool TriLinearHexahedronTransformation< DataType >::inverse(
    DataType x_phy, DataType y_phy, DataType z_phy, DataType &x_ref,
    DataType &y_ref, DataType &z_ref) const {
  bool success =
      this->inverse_by_decomposition(x_phy, y_phy, z_phy, x_ref, y_ref, z_ref);

  bool found = false;
  if (success) {
    // reference point obtained by decomposition is sufficiently accurate ->
    // exit
    return true;
  } else {
    // reference point obtained by decomposition is not accurate enough -> use
    // as initial value for Newton method
    found = this->inverse_newton_3d(x_phy, y_phy, z_phy, x_ref, y_ref, z_ref,
                                    x_ref, y_ref, z_ref);
  }

#ifndef NDEBUG
  Coord coord_ref(3);
  coord_ref[0] = x_ref;
  coord_ref[1] = y_ref;
  coord_ref[2] = z_ref;

  DataType x = this->x(coord_ref);
  DataType y = this->y(coord_ref);
  DataType z = this->z(coord_ref);

  // assert ( std::abs ( x - x_phy ) < 1e-12 );
  // assert ( std::abs ( y - y_phy ) < 1e-12 );
  // assert ( std::abs ( z - z_phy ) < 1e-12 );
  /*
     std::cout << " a inverse " << std::endl;

     if ( !found )
     {
         std::cout << " algined inverse: " << x_phy << ", " << y_phy << ", " <<
     z_phy << " -> " << x_ref << ", " << y_ref << ", " << z_ref << std::endl;
         std::cout << this->coord_vtx_[this->ij2ind ( 0, 0 )] << ", " <<
     this->coord_vtx_[this->ij2ind ( 0, 1 )] << ", " <<
     this->coord_vtx_[this->ij2ind ( 0, 2 )] << std::endl; std::cout <<
     this->coord_vtx_[this->ij2ind ( 1, 0 )] << ", " <<
     this->coord_vtx_[this->ij2ind ( 1, 1 )] << ", " <<
     this->coord_vtx_[this->ij2ind ( 1, 2 )] << std::endl; std::cout <<
     this->coord_vtx_[this->ij2ind ( 2, 0 )] << ", " <<
     this->coord_vtx_[this->ij2ind ( 2, 1 )] << ", " <<
     this->coord_vtx_[this->ij2ind ( 2, 2 )] << std::endl; std::cout <<
     this->coord_vtx_[this->ij2ind ( 3, 0 )] << ", " <<
     this->coord_vtx_[this->ij2ind ( 3, 1 )] << ", " <<
     this->coord_vtx_[this->ij2ind ( 3, 2 )] << std::endl; std::cout <<
     this->coord_vtx_[this->ij2ind ( 4, 0 )] << ", " <<
     this->coord_vtx_[this->ij2ind ( 4, 1 )] << ", " <<
     this->coord_vtx_[this->ij2ind ( 4, 2 )] << std::endl; std::cout <<
     this->coord_vtx_[this->ij2ind ( 5, 0 )] << ", " <<
     this->coord_vtx_[this->ij2ind ( 5, 1 )] << ", " <<
     this->coord_vtx_[this->ij2ind ( 5, 2 )] << std::endl; std::cout <<
     this->coord_vtx_[this->ij2ind ( 6, 0 )] << ", " <<
     this->coord_vtx_[this->ij2ind ( 6, 1 )] << ", " <<
     this->coord_vtx_[this->ij2ind ( 6, 2 )] << std::endl; std::cout <<
     this->coord_vtx_[this->ij2ind ( 7, 0 )] << ", " <<
     this->coord_vtx_[this->ij2ind ( 7, 1 )] << ", " <<
     this->coord_vtx_[this->ij2ind ( 7, 2 )] << std::endl;
     }*/
#endif
  return found;

  // Newton method with old initial value
  // return this->inverse_newton_3d ( x_phy, y_phy, z_phy, x_ref, y_ref, z_ref
  // );
}

template < class DataType >
DataType
TriLinearHexahedronTransformation< DataType >::x(const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_1 = coord_ref[1];
  const DataType coord_2 = coord_ref[2];

  return +this->coord_vtx_[this->ij2ind(0, 0)] * (1. - coord_0) *
             (1. - coord_1) * (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(1, 0)] * coord_0 * (1. - coord_1) *
             (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(2, 0)] * coord_0 * coord_1 *
             (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(3, 0)] * (1. - coord_0) * coord_1 *
             (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(4, 0)] * (1. - coord_0) *
             (1. - coord_1) * coord_2 +
         this->coord_vtx_[this->ij2ind(5, 0)] * coord_0 * (1. - coord_1) *
             coord_2 +
         this->coord_vtx_[this->ij2ind(6, 0)] * coord_0 * coord_1 * coord_2 +
         this->coord_vtx_[this->ij2ind(7, 0)] * (1. - coord_0) * coord_1 *
             coord_2;
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::x_x(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  const DataType coord_1 = coord_ref[1];
  const DataType coord_2 = coord_ref[2];

  return -this->coord_vtx_[this->ij2ind(0, 0)] * (1. - coord_1) *
             (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(1, 0)] * (1. - coord_1) *
             (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(2, 0)] * coord_1 * (1. - coord_2) -
         this->coord_vtx_[this->ij2ind(3, 0)] * coord_1 * (1. - coord_2) -
         this->coord_vtx_[this->ij2ind(4, 0)] * (1. - coord_1) * coord_2 +
         this->coord_vtx_[this->ij2ind(5, 0)] * (1. - coord_1) * coord_2 +
         this->coord_vtx_[this->ij2ind(6, 0)] * coord_1 * coord_2 -
         this->coord_vtx_[this->ij2ind(7, 0)] * coord_1 * coord_2;
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::x_y(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_2 = coord_ref[2];

  return -this->coord_vtx_[this->ij2ind(0, 0)] * (1. - coord_0) *
             (1. - coord_2) -
         this->coord_vtx_[this->ij2ind(1, 0)] * coord_0 * (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(2, 0)] * coord_0 * (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(3, 0)] * (1. - coord_0) *
             (1. - coord_2) -
         this->coord_vtx_[this->ij2ind(4, 0)] * (1. - coord_0) * coord_2 -
         this->coord_vtx_[this->ij2ind(5, 0)] * coord_0 * coord_2 +
         this->coord_vtx_[this->ij2ind(6, 0)] * coord_0 * coord_2 +
         this->coord_vtx_[this->ij2ind(7, 0)] * (1. - coord_0) * coord_2;
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::x_z(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_1 = coord_ref[1];

  return -this->coord_vtx_[this->ij2ind(0, 0)] * (1. - coord_0) *
             (1. - coord_1) -
         this->coord_vtx_[this->ij2ind(1, 0)] * coord_0 * (1. - coord_1) -
         this->coord_vtx_[this->ij2ind(2, 0)] * coord_0 * coord_1 -
         this->coord_vtx_[this->ij2ind(3, 0)] * (1. - coord_0) * coord_1 +
         this->coord_vtx_[this->ij2ind(4, 0)] * (1. - coord_0) *
             (1. - coord_1) +
         this->coord_vtx_[this->ij2ind(5, 0)] * coord_0 * (1. - coord_1) +
         this->coord_vtx_[this->ij2ind(6, 0)] * coord_0 * coord_1 +
         this->coord_vtx_[this->ij2ind(7, 0)] * (1. - coord_0) * coord_1;
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::x_xx(
    const Coord &coord_ref) const {
  return 0.;
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::x_xy(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  const DataType coord_2 = coord_ref[2];

  return +this->coord_vtx_[this->ij2ind(0, 0)] * (1. - coord_2) -
         this->coord_vtx_[this->ij2ind(1, 0)] * (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(2, 0)] * (1. - coord_2) -
         this->coord_vtx_[this->ij2ind(3, 0)] * (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(4, 0)] * coord_2 -
         this->coord_vtx_[this->ij2ind(5, 0)] * coord_2 +
         this->coord_vtx_[this->ij2ind(6, 0)] * coord_2 -
         this->coord_vtx_[this->ij2ind(7, 0)] * coord_2;
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::x_xz(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  const DataType coord_1 = coord_ref[1];

  return +this->coord_vtx_[this->ij2ind(0, 0)] * (1. - coord_1) -
         this->coord_vtx_[this->ij2ind(1, 0)] * (1. - coord_1) -
         this->coord_vtx_[this->ij2ind(2, 0)] * coord_1 +
         this->coord_vtx_[this->ij2ind(3, 0)] * coord_1 -
         this->coord_vtx_[this->ij2ind(4, 0)] * (1. - coord_1) *
             +this->coord_vtx_[this->ij2ind(5, 0)] * (1. - coord_1) *
             +this->coord_vtx_[this->ij2ind(6, 0)] * coord_1 -
         this->coord_vtx_[this->ij2ind(7, 0)] * coord_1;
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::x_yy(
    const Coord &coord_ref) const {
  return 0.;
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::x_yz(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  const DataType coord_0 = coord_ref[0];

  return +this->coord_vtx_[this->ij2ind(0, 0)] * (1. - coord_0) +
         this->coord_vtx_[this->ij2ind(1, 0)] * coord_0 -
         this->coord_vtx_[this->ij2ind(2, 0)] * coord_0 -
         this->coord_vtx_[this->ij2ind(3, 0)] * (1. - coord_0) -
         this->coord_vtx_[this->ij2ind(4, 0)] * (1. - coord_0) -
         this->coord_vtx_[this->ij2ind(5, 0)] * coord_0 +
         this->coord_vtx_[this->ij2ind(6, 0)] * coord_0 +
         this->coord_vtx_[this->ij2ind(7, 0)] * (1. - coord_0);
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::x_zz(
    const Coord &coord_ref) const {
  return 0.;
}

template < class DataType >
DataType
TriLinearHexahedronTransformation< DataType >::y(const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_1 = coord_ref[1];
  const DataType coord_2 = coord_ref[2];

  return +this->coord_vtx_[this->ij2ind(0, 1)] * (1. - coord_0) *
             (1. - coord_1) * (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(1, 1)] * coord_0 * (1. - coord_1) *
             (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(2, 1)] * coord_0 * coord_1 *
             (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(3, 1)] * (1. - coord_0) * coord_1 *
             (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(4, 1)] * (1. - coord_0) *
             (1. - coord_1) * coord_2 +
         this->coord_vtx_[this->ij2ind(5, 1)] * coord_0 * (1. - coord_1) *
             coord_2 +
         this->coord_vtx_[this->ij2ind(6, 1)] * coord_0 * coord_1 * coord_2 +
         this->coord_vtx_[this->ij2ind(7, 1)] * (1. - coord_0) * coord_1 *
             coord_2;
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::y_x(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  const DataType coord_1 = coord_ref[1];
  const DataType coord_2 = coord_ref[2];

  return -this->coord_vtx_[this->ij2ind(0, 1)] * (1. - coord_1) *
             (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(1, 1)] * (1. - coord_1) *
             (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(2, 1)] * coord_1 * (1. - coord_2) -
         this->coord_vtx_[this->ij2ind(3, 1)] * coord_1 * (1. - coord_2) -
         this->coord_vtx_[this->ij2ind(4, 1)] * (1. - coord_1) * coord_2 +
         this->coord_vtx_[this->ij2ind(5, 1)] * (1. - coord_1) * coord_2 +
         this->coord_vtx_[this->ij2ind(6, 1)] * coord_1 * coord_2 -
         this->coord_vtx_[this->ij2ind(7, 1)] * coord_1 * coord_2;
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::y_y(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_2 = coord_ref[2];

  return -this->coord_vtx_[this->ij2ind(0, 1)] * (1. - coord_0) *
             (1. - coord_2) -
         this->coord_vtx_[this->ij2ind(1, 1)] * coord_0 * (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(2, 1)] * coord_0 * (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(3, 1)] * (1. - coord_0) *
             (1. - coord_2) -
         this->coord_vtx_[this->ij2ind(4, 1)] * (1. - coord_0) * coord_2 -
         this->coord_vtx_[this->ij2ind(5, 1)] * coord_0 * coord_2 +
         this->coord_vtx_[this->ij2ind(6, 1)] * coord_0 * coord_2 +
         this->coord_vtx_[this->ij2ind(7, 1)] * (1. - coord_0) * coord_2;
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::y_z(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_1 = coord_ref[1];

  return -this->coord_vtx_[this->ij2ind(0, 1)] * (1. - coord_0) *
             (1. - coord_1) -
         this->coord_vtx_[this->ij2ind(1, 1)] * coord_0 * (1. - coord_1) -
         this->coord_vtx_[this->ij2ind(2, 1)] * coord_0 * coord_1 -
         this->coord_vtx_[this->ij2ind(3, 1)] * (1. - coord_0) * coord_1 +
         this->coord_vtx_[this->ij2ind(4, 1)] * (1. - coord_0) *
             (1. - coord_1) +
         this->coord_vtx_[this->ij2ind(5, 1)] * coord_0 * (1. - coord_1) +
         this->coord_vtx_[this->ij2ind(6, 1)] * coord_0 * coord_1 +
         this->coord_vtx_[this->ij2ind(7, 1)] * (1. - coord_0) * coord_1;
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::y_xx(
    const Coord &coord_ref) const {
  return 0.;
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::y_xy(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  const DataType coord_2 = coord_ref[2];

  return +this->coord_vtx_[this->ij2ind(0, 1)] * (1. - coord_2) -
         this->coord_vtx_[this->ij2ind(1, 1)] * (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(2, 1)] * (1. - coord_2) -
         this->coord_vtx_[this->ij2ind(3, 1)] * (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(4, 1)] * coord_2 -
         this->coord_vtx_[this->ij2ind(5, 1)] * coord_2 +
         this->coord_vtx_[this->ij2ind(6, 1)] * coord_2 -
         this->coord_vtx_[this->ij2ind(7, 1)] * coord_2;
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::y_xz(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  const DataType coord_1 = coord_ref[1];

  return +this->coord_vtx_[this->ij2ind(0, 1)] * (1. - coord_1) -
         this->coord_vtx_[this->ij2ind(1, 1)] * (1. - coord_1) -
         this->coord_vtx_[this->ij2ind(2, 1)] * coord_1 +
         this->coord_vtx_[this->ij2ind(3, 1)] * coord_1 -
         this->coord_vtx_[this->ij2ind(4, 1)] * (1. - coord_1) +
         this->coord_vtx_[this->ij2ind(5, 1)] * (1. - coord_1) +
         this->coord_vtx_[this->ij2ind(6, 1)] * coord_1 -
         this->coord_vtx_[this->ij2ind(7, 1)] * coord_1;
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::y_yy(
    const Coord &coord_ref) const {
  return 0.;
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::y_yz(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  const DataType coord_0 = coord_ref[0];

  return +this->coord_vtx_[this->ij2ind(0, 1)] * (1. - coord_0) +
         this->coord_vtx_[this->ij2ind(1, 1)] * coord_0 -
         this->coord_vtx_[this->ij2ind(2, 1)] * coord_0 -
         this->coord_vtx_[this->ij2ind(3, 1)] * (1. - coord_0) -
         this->coord_vtx_[this->ij2ind(4, 1)] * (1. - coord_0) -
         this->coord_vtx_[this->ij2ind(5, 1)] * coord_0 +
         this->coord_vtx_[this->ij2ind(6, 1)] * coord_0 +
         this->coord_vtx_[this->ij2ind(7, 1)] * (1. - coord_0);
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::y_zz(
    const Coord &coord_ref) const {
  return 0.;
}

template < class DataType >
DataType
TriLinearHexahedronTransformation< DataType >::z(const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_1 = coord_ref[1];
  const DataType coord_2 = coord_ref[2];

  return +this->coord_vtx_[this->ij2ind(0, 2)] * (1. - coord_0) *
             (1. - coord_1) * (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(1, 2)] * coord_0 * (1. - coord_1) *
             (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(2, 2)] * coord_0 * coord_1 *
             (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(3, 2)] * (1. - coord_0) * coord_1 *
             (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(4, 2)] * (1. - coord_0) *
             (1. - coord_1) * coord_2 +
         this->coord_vtx_[this->ij2ind(5, 2)] * coord_0 * (1. - coord_1) *
             coord_2 +
         this->coord_vtx_[this->ij2ind(6, 2)] * coord_0 * coord_1 * coord_2 +
         this->coord_vtx_[this->ij2ind(7, 2)] * (1. - coord_0) * coord_1 *
             coord_2;
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::z_x(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  const DataType coord_1 = coord_ref[1];
  const DataType coord_2 = coord_ref[2];

  return -this->coord_vtx_[this->ij2ind(0, 2)] * (1. - coord_1) *
             (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(1, 2)] * (1. - coord_1) *
             (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(2, 2)] * coord_1 * (1. - coord_2) -
         this->coord_vtx_[this->ij2ind(3, 2)] * coord_1 * (1. - coord_2) -
         this->coord_vtx_[this->ij2ind(4, 2)] * (1. - coord_1) * coord_2 +
         this->coord_vtx_[this->ij2ind(5, 2)] * (1. - coord_1) * coord_2 +
         this->coord_vtx_[this->ij2ind(6, 2)] * coord_1 * coord_2 -
         this->coord_vtx_[this->ij2ind(7, 2)] * coord_1 * coord_2;
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::z_y(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_2 = coord_ref[2];

  return -this->coord_vtx_[this->ij2ind(0, 2)] * (1. - coord_0) *
             (1. - coord_2) -
         this->coord_vtx_[this->ij2ind(1, 2)] * coord_0 * (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(2, 2)] * coord_0 * (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(3, 2)] * (1. - coord_0) *
             (1. - coord_2) -
         this->coord_vtx_[this->ij2ind(4, 2)] * (1. - coord_0) * coord_2 -
         this->coord_vtx_[this->ij2ind(5, 2)] * coord_0 * coord_2 +
         this->coord_vtx_[this->ij2ind(6, 2)] * coord_0 * coord_2 +
         this->coord_vtx_[this->ij2ind(7, 2)] * (1. - coord_0) * coord_2;
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::z_z(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_1 = coord_ref[1];

  return -this->coord_vtx_[this->ij2ind(0, 2)] * (1. - coord_0) *
             (1. - coord_1) -
         this->coord_vtx_[this->ij2ind(1, 2)] * coord_0 * (1. - coord_1) -
         this->coord_vtx_[this->ij2ind(2, 2)] * coord_0 * coord_1 -
         this->coord_vtx_[this->ij2ind(3, 2)] * (1. - coord_0) * coord_1 +
         this->coord_vtx_[this->ij2ind(4, 2)] * (1. - coord_0) *
             (1. - coord_1) +
         this->coord_vtx_[this->ij2ind(5, 2)] * coord_0 * (1. - coord_1) +
         this->coord_vtx_[this->ij2ind(6, 2)] * coord_0 * coord_1 +
         this->coord_vtx_[this->ij2ind(7, 2)] * (1. - coord_0) * coord_1;
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::z_xx(
    const Coord &coord_ref) const {
  return 0.;
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::z_xy(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  const DataType coord_2 = coord_ref[2];

  return +this->coord_vtx_[this->ij2ind(0, 2)] * (1. - coord_2) -
         this->coord_vtx_[this->ij2ind(1, 2)] * (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(2, 2)] * (1. - coord_2) -
         this->coord_vtx_[this->ij2ind(3, 2)] * (1. - coord_2) +
         this->coord_vtx_[this->ij2ind(4, 2)] * coord_2 -
         this->coord_vtx_[this->ij2ind(5, 2)] * coord_2 +
         this->coord_vtx_[this->ij2ind(6, 2)] * coord_2 -
         this->coord_vtx_[this->ij2ind(7, 2)] * coord_2;
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::z_xz(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  const DataType coord_1 = coord_ref[1];

  return +this->coord_vtx_[this->ij2ind(0, 2)] * (1. - coord_1) -
         this->coord_vtx_[this->ij2ind(1, 2)] * (1. - coord_1) -
         this->coord_vtx_[this->ij2ind(2, 2)] * coord_1 +
         this->coord_vtx_[this->ij2ind(3, 2)] * coord_1 -
         this->coord_vtx_[this->ij2ind(4, 2)] * (1. - coord_1) +
         this->coord_vtx_[this->ij2ind(5, 2)] * (1. - coord_1) +
         this->coord_vtx_[this->ij2ind(6, 2)] * coord_1 -
         this->coord_vtx_[this->ij2ind(7, 2)] * coord_1;
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::z_yy(
    const Coord &coord_ref) const {
  return 0.;
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::z_yz(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == 3);
  const DataType coord_0 = coord_ref[0];

  return +this->coord_vtx_[this->ij2ind(0, 2)] * (1. - coord_0) +
         this->coord_vtx_[this->ij2ind(1, 2)] * coord_0 -
         this->coord_vtx_[this->ij2ind(2, 2)] * coord_0 -
         this->coord_vtx_[this->ij2ind(3, 2)] * (1. - coord_0) -
         this->coord_vtx_[this->ij2ind(4, 2)] * (1. - coord_0) -
         this->coord_vtx_[this->ij2ind(5, 2)] * coord_0 +
         this->coord_vtx_[this->ij2ind(6, 2)] * coord_0 +
         this->coord_vtx_[this->ij2ind(7, 2)] * (1. - coord_0);
}

template < class DataType >
DataType TriLinearHexahedronTransformation< DataType >::z_zz(
    const Coord &coord_ref) const {
  return 0.;
}

template < class DataType >
bool TriLinearHexahedronTransformation< DataType >::contains_reference_point(
    const Coord &coord_ref) const {
  assert(coord_ref.size() == this->gdim_);
  return coord_ref[0] >= -this->eps_ && coord_ref[0] <= 1. + this->eps_ &&
         coord_ref[1] >= -this->eps_ && coord_ref[1] <= 1. + this->eps_ &&
         coord_ref[2] >= -this->eps_ && coord_ref[2] <= 1. + this->eps_;
}

template class TriLinearHexahedronTransformation< double >;
template class TriLinearHexahedronTransformation< float >;

} // namespace doffem
} // namespace hiflow
