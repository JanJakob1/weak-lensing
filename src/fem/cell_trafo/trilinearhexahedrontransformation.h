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

#ifndef __FEM_TRILINEAR_HEXAHEDRON_TRANSFORMATION_H_
#define __FEM_TRILINEAR_HEXAHEDRON_TRANSFORMATION_H_

#include "common/log.h"
#include "fem/cell_trafo/cell_transformation.h"
#include "fem/cell_trafo/lineartetrahedrontransformation.h"
#include <cmath>
#include <iomanip>

namespace hiflow {
namespace doffem {

///
/// \class TriLinearHexahedronTransformation trilinearhexahedrontransformation.h
/// \brief Trilinear transformation mapping from reference to physical cell for
/// a Hexahedron \author \author Michael Schick<br>Martin Baumann<br>Simon Gawlok<br>Philipp Gerstner
///

template < class DataType, int DIM >
class TriLinearHexahedronTransformation
    : public CellTransformation< DataType, DIM > {
public:
  typedef Vec<DIM, DataType> Coord;

  explicit TriLinearHexahedronTransformation(ConstRefCellPtr<DataType, DIM> ref_cell);

  bool differs_by_translation_from (ConstCellTrafoPtr<DataType, DIM> rhs) const;
   
  bool inverse(Coord co_phy, Coord &co_ref) const;

  DataType x(const Coord &coord_ref) const;
  DataType x_x(const Coord &coord_ref) const;
  DataType x_y(const Coord &coord_ref) const;
  DataType x_z(const Coord &coord_ref) const;
  DataType x_xy(const Coord &coord_ref) const;
  DataType x_xz(const Coord &coord_ref) const;
  DataType x_yz(const Coord &coord_ref) const;
  DataType x_xyz(const Coord &coord_ref) const;
  
  DataType y(const Coord &coord_ref) const;
  DataType y_x(const Coord &coord_ref) const;
  DataType y_y(const Coord &coord_ref) const;
  DataType y_z(const Coord &coord_ref) const;
  DataType y_xy(const Coord &coord_ref) const;
  DataType y_xz(const Coord &coord_ref) const;
  DataType y_yz(const Coord &coord_ref) const;
  DataType y_xyz(const Coord &coord_ref) const;
  
  DataType z(const Coord &coord_ref) const;
  DataType z_x(const Coord &coord_ref) const;
  DataType z_y(const Coord &coord_ref) const;
  DataType z_z(const Coord &coord_ref) const;
  DataType z_xy(const Coord &coord_ref) const;
  DataType z_xz(const Coord &coord_ref) const;
  DataType z_yz(const Coord &coord_ref) const;
  DataType z_xyz(const Coord &coord_ref) const;
  
protected:
  bool inverse_by_decomposition(Coord co_phy, Coord &co_ref) const;

  int8_t tetrahedron_decomp_ind_[5][4];
};

// Reordering of vertices to make transformation coorespond to mesh
// ordering, with (0,0,0) mapped to vertex 0, and (1,1,1) mapped to vertex 7.

template < class DataType, int DIM >
TriLinearHexahedronTransformation<DataType, DIM >::TriLinearHexahedronTransformation(ConstRefCellPtr<DataType, DIM> ref_cell)
    : CellTransformation< DataType, DIM >(ref_cell) 
{
  this->order_ = 3;  
  this->fixed_ref_cell_type_ = REF_CELL_HEX_STD;
  this->name_ = "Hex";
   

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

template < class DataType, int DIM >
bool TriLinearHexahedronTransformation< DataType, DIM >::differs_by_translation_from (ConstCellTrafoPtr<DataType, DIM> rhs) const
{
  assert(this->ref_cell_);
  
  if (this->ref_cell_->type() != rhs->get_ref_cell()->type())
  {
    return false;
  }
  
  std::vector< Coord > ref_coords = this->ref_cell_->get_coords();
  assert (ref_coords.size() == 8);
  
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
  Coord my_p45 = this->transform(ref_coords[4]) - this->transform(ref_coords[5]); 
  Coord rhs_p45 = rhs->transform(ref_coords[4]) - rhs->transform(ref_coords[5]); 
  if (my_p45 != rhs_p45)
  {
    return false;
  }
  Coord my_p56 = this->transform(ref_coords[5]) - this->transform(ref_coords[6]); 
  Coord rhs_p56 = rhs->transform(ref_coords[5]) - rhs->transform(ref_coords[6]); 
  if (my_p56 != rhs_p56)
  {
    return false;
  }
  Coord my_p67 = this->transform(ref_coords[6]) - this->transform(ref_coords[7]); 
  Coord rhs_p67 = rhs->transform(ref_coords[6]) - rhs->transform(ref_coords[7]); 
  if (my_p67 != rhs_p67)
  {
    return false;
  }
  Coord my_p74 = this->transform(ref_coords[7]) - this->transform(ref_coords[4]); 
  Coord rhs_p74 = rhs->transform(ref_coords[7]) - rhs->transform(ref_coords[4]); 
  if (my_p74 != rhs_p74)
  {
    return false;
  }
  Coord my_p04 = this->transform(ref_coords[0]) - this->transform(ref_coords[4]); 
  Coord rhs_p04 = rhs->transform(ref_coords[0]) - rhs->transform(ref_coords[4]); 
  if (my_p04 != rhs_p04)
  {
    return false;
  }
  Coord my_p15 = this->transform(ref_coords[1]) - this->transform(ref_coords[5]); 
  Coord rhs_p15 = rhs->transform(ref_coords[1]) - rhs->transform(ref_coords[5]); 
  if (my_p15 != rhs_p15)
  {
    return false;
  }
  Coord my_p26 = this->transform(ref_coords[2]) - this->transform(ref_coords[6]); 
  Coord rhs_p26 = rhs->transform(ref_coords[2]) - rhs->transform(ref_coords[6]); 
  if (my_p26 != rhs_p26)
  {
    return false;
  }
  Coord my_p37 = this->transform(ref_coords[3]) - this->transform(ref_coords[7]); 
  Coord rhs_p37 = rhs->transform(ref_coords[3]) - rhs->transform(ref_coords[7]); 
  if (my_p37 != rhs_p37)
  {
    return false;
  }
  return true;
}

template < class DataType, int DIM >
bool TriLinearHexahedronTransformation< DataType, DIM >::inverse_by_decomposition( Coord co_phy, Coord &co_ref) const 
{
  assert (DIM == 3);
  RefCellPtr<DataType, DIM> ref_cell_tet = RefCellPtr<DataType, DIM>(new RefCellTetStd<DataType, DIM> );
  LinearTetrahedronTransformation< DataType, DIM > tetra_trafo(ref_cell_tet);
  Coord ref_tet;

  int tetra_ind = -1;
  
  // loop through tetrahedron decomposition
  for (int t = 0; t < 5; ++t) 
  {
    // build linear tetrahedron transformation
    std::vector<DataType> tetra_coord(12);

    for (int v = 0; v < 4; ++v) 
    {
      for (int d = 0; d < 3; ++d) 
      {
        tetra_coord[this->ij2ind(v,d)] = this->coord_vtx_[this->tetrahedron_decomp_ind_[t][v]][d];
      }
    }
    tetra_trafo.reinit(tetra_coord);

    // compute reference coordinates w.r.t. tetrahedron
    bool tetra_success = tetra_trafo.inverse(co_phy, ref_tet);
    if (!tetra_success) 
    {
      continue;
    }

    // check whether reference coordinates are contained in tetrahedron
    // if yes: convert reference coordinates
    if (tetra_trafo.contains_reference_point(ref_tet)) 
    {
      switch (t) 
      {
      case 0:
        co_ref[0] = 1. - ref_tet[0];
        co_ref[1] = 1. - ref_tet[1];
        co_ref[2] = ref_tet[2];
        break;
      case 1:
        co_ref[0] = ref_tet[0];
        co_ref[1] = ref_tet[1];
        co_ref[2] = ref_tet[2];
        break;
      case 2:
        co_ref[0] = 1. - ref_tet[0];
        co_ref[1] = ref_tet[1];
        co_ref[2] = 1. - ref_tet[2];
        break;
      case 3:
        co_ref[0] = ref_tet[0];
        co_ref[1] = 1. - ref_tet[1];
        co_ref[2] = 1. - ref_tet[2];
        break;
      case 4:
        co_ref[0] = 1. - ref_tet[0] - ref_tet[2];
        co_ref[1] = 1. - ref_tet[1] - ref_tet[2];
        co_ref[2] = 1. - ref_tet[0] - ref_tet[1];
        break;
      }
      tetra_ind = t;
      break;
    }
  }
  // TODO: is this check necessary?
  // check if reference point found by decomposition is sufficiently accurate
  Coord phy_test;
  phy_test[0] = this->x(co_ref);
  phy_test[1] = this->y(co_ref);
  phy_test[2] = this->z(co_ref);

  DataType residual = norm (phy_test - co_phy);

#ifndef NDEBUG
  LOG_DEBUG(3, "==============");
  LOG_DEBUG(3, "Physical point " << co_phy[0] << " " << co_phy[1] << " " << co_phy[2]
                                 << " is contained in tetrahedron "
                                 << tetra_ind);

  if (tetra_ind == -1) {
    LOG_DEBUG(3, "Physical point " << co_phy[0] << " " << co_phy[1] << " " << co_phy[2]
                                   << " is not contained in any tetrahedron ");
    for (int v = 0; v < 8; ++v) {
      LOG_DEBUG(3, "vertex " << v << " : "
                             << this->coord_vtx_[v][0] << ", "
                             << this->coord_vtx_[v][1] << ", "
                             << this->coord_vtx_[v][2]);
    }
  }

  LOG_DEBUG(3, "Inverse cell-trafo based on tetrahedron decomposition yields residual = " << residual);
#endif

  if (residual < 5. * std::numeric_limits< DataType >::epsilon()) {
    return true;
  } 
  else 
  {
    return false;
  }
}

template < class DataType, int DIM >
bool TriLinearHexahedronTransformation< DataType, DIM >::inverse(Coord co_phy, Coord &co_ref) const 
{
  assert (DIM == 3);
  bool success = this->inverse_by_decomposition(co_phy, co_ref);

  bool found = false;
  if (success) 
  {
    // reference point obtained by decomposition is sufficiently accurate -> exit
    return true;
  } 
  else 
  {
    // reference point obtained by decomposition is not accurate enough -> use
    // as initial value for Newton method
    bool newton_success = this->inverse_newton(co_phy, co_ref, co_ref);
    if (newton_success)
    {
      return this->contains_reference_point(co_ref);
    }
    else
    {
      return false;
    }
  }

#if 0
  DataType x = this->x(co_ref);
  DataType y = this->y(co_ref);
  DataType z = this->z(co_ref);

  assert ( std::abs ( x - x_phy ) < 1e-12 );
  assert ( std::abs ( y - y_phy ) < 1e-12 );
  assert ( std::abs ( z - z_phy ) < 1e-12 );
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
  // Newton method with old initial value
  // return this->inverse_newton_3d ( x_phy, y_phy, z_phy, x_ref, y_ref, z_ref
  // );
}

template < class DataType, int DIM >
DataType
TriLinearHexahedronTransformation< DataType, DIM >::x(const Coord &coord_ref) const {
  assert (DIM == 3);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_1 = coord_ref[1];
  const DataType coord_2 = coord_ref[2];

  return +this->coord_vtx_[0][0] * (1. - coord_0) *
             (1. - coord_1) * (1. - coord_2) +
         this->coord_vtx_[1][0] * coord_0 * (1. - coord_1) *
             (1. - coord_2) +
         this->coord_vtx_[2][0] * coord_0 * coord_1 *
             (1. - coord_2) +
         this->coord_vtx_[3][0] * (1. - coord_0) * coord_1 *
             (1. - coord_2) +
         this->coord_vtx_[4][0] * (1. - coord_0) *
             (1. - coord_1) * coord_2 +
         this->coord_vtx_[5][0] * coord_0 * (1. - coord_1) *
             coord_2 +
         this->coord_vtx_[6][0] * coord_0 * coord_1 * coord_2 +
         this->coord_vtx_[7][0] * (1. - coord_0) * coord_1 *
             coord_2;
}

template < class DataType, int DIM >
DataType TriLinearHexahedronTransformation< DataType, DIM >::x_x( const Coord &coord_ref) const {
  assert (DIM == 3);
  const DataType coord_1 = coord_ref[1];
  const DataType coord_2 = coord_ref[2];

  return -this->coord_vtx_[0][0] * (1. - coord_1) *
             (1. - coord_2) +
         this->coord_vtx_[1][0] * (1. - coord_1) *
             (1. - coord_2) +
         this->coord_vtx_[2][0] * coord_1 * (1. - coord_2) -
         this->coord_vtx_[3][0] * coord_1 * (1. - coord_2) -
         this->coord_vtx_[4][0] * (1. - coord_1) * coord_2 +
         this->coord_vtx_[5][0] * (1. - coord_1) * coord_2 +
         this->coord_vtx_[6][0] * coord_1 * coord_2 -
         this->coord_vtx_[7][0] * coord_1 * coord_2;
}

template < class DataType, int DIM >
DataType TriLinearHexahedronTransformation< DataType, DIM >::x_y( const Coord &coord_ref) const {
  assert (DIM == 3);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_2 = coord_ref[2];

  return -this->coord_vtx_[0][0] * (1. - coord_0) *
             (1. - coord_2) -
         this->coord_vtx_[1][0] * coord_0 * (1. - coord_2) +
         this->coord_vtx_[2][0] * coord_0 * (1. - coord_2) +
         this->coord_vtx_[3][0] * (1. - coord_0) *
             (1. - coord_2) -
         this->coord_vtx_[4][0] * (1. - coord_0) * coord_2 -
         this->coord_vtx_[5][0] * coord_0 * coord_2 +
         this->coord_vtx_[6][0] * coord_0 * coord_2 +
         this->coord_vtx_[7][0] * (1. - coord_0) * coord_2;
}

template < class DataType, int DIM >
DataType TriLinearHexahedronTransformation< DataType, DIM >::x_z( const Coord &coord_ref) const {
  assert (DIM == 3);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_1 = coord_ref[1];

  return -this->coord_vtx_[0][0] * (1. - coord_0) *
             (1. - coord_1) -
         this->coord_vtx_[1][0] * coord_0 * (1. - coord_1) -
         this->coord_vtx_[2][0] * coord_0 * coord_1 -
         this->coord_vtx_[3][0] * (1. - coord_0) * coord_1 +
         this->coord_vtx_[4][0] * (1. - coord_0) *
             (1. - coord_1) +
         this->coord_vtx_[5][0] * coord_0 * (1. - coord_1) +
         this->coord_vtx_[6][0] * coord_0 * coord_1 +
         this->coord_vtx_[7][0] * (1. - coord_0) * coord_1;
}

template < class DataType, int DIM >
DataType TriLinearHexahedronTransformation< DataType, DIM >::x_xy( const Coord &coord_ref) const {
  assert (DIM == 3);
  const DataType coord_2 = coord_ref[2];

  return +this->coord_vtx_[0][0] * (1. - coord_2) -
         this->coord_vtx_[1][0] * (1. - coord_2) +
         this->coord_vtx_[2][0] * (1. - coord_2) -
         this->coord_vtx_[3][0] * (1. - coord_2) +
         this->coord_vtx_[4][0] * coord_2 -
         this->coord_vtx_[5][0] * coord_2 +
         this->coord_vtx_[6][0] * coord_2 -
         this->coord_vtx_[7][0] * coord_2;
}

template < class DataType, int DIM >
DataType TriLinearHexahedronTransformation< DataType, DIM >::x_xz( const Coord &coord_ref) const {
  assert (DIM == 3);
  const DataType coord_1 = coord_ref[1];

  return +this->coord_vtx_[0][0] * (1. - coord_1) -
         this->coord_vtx_[1][0] * (1. - coord_1) -
         this->coord_vtx_[2][0] * coord_1 +
         this->coord_vtx_[3][0] * coord_1 -
         this->coord_vtx_[4][0] * (1. - coord_1) *
             +this->coord_vtx_[5][0] * (1. - coord_1) *
             +this->coord_vtx_[6][0] * coord_1 -
         this->coord_vtx_[7][0] * coord_1;
}

template < class DataType, int DIM >
DataType TriLinearHexahedronTransformation< DataType, DIM >::x_yz( const Coord &coord_ref) const {
  assert (DIM == 3);
  const DataType coord_0 = coord_ref[0];

  return +this->coord_vtx_[0][0] * (1. - coord_0) +
         this->coord_vtx_[1][0] * coord_0 -
         this->coord_vtx_[2][0] * coord_0 -
         this->coord_vtx_[3][0] * (1. - coord_0) -
         this->coord_vtx_[4][0] * (1. - coord_0) -
         this->coord_vtx_[5][0] * coord_0 +
         this->coord_vtx_[6][0] * coord_0 +
         this->coord_vtx_[7][0] * (1. - coord_0);
}

template < class DataType, int DIM >
DataType TriLinearHexahedronTransformation< DataType, DIM >::x_xyz( const Coord &coord_ref) const {
  assert (DIM == 3);

  return - this->coord_vtx_[0][0]
         + this->coord_vtx_[1][0] 
         - this->coord_vtx_[2][0] 
         + this->coord_vtx_[3][0] 
         + this->coord_vtx_[4][0] 
         - this->coord_vtx_[5][0] 
         + this->coord_vtx_[6][0]
         - this->coord_vtx_[7][0];
}

template < class DataType, int DIM >
DataType
TriLinearHexahedronTransformation< DataType, DIM >::y(const Coord &coord_ref) const {
  assert (DIM == 3);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_1 = coord_ref[1];
  const DataType coord_2 = coord_ref[2];

  return +this->coord_vtx_[0][1] * (1. - coord_0) *
             (1. - coord_1) * (1. - coord_2) +
         this->coord_vtx_[1][1] * coord_0 * (1. - coord_1) *
             (1. - coord_2) +
         this->coord_vtx_[2][1] * coord_0 * coord_1 *
             (1. - coord_2) +
         this->coord_vtx_[3][1] * (1. - coord_0) * coord_1 *
             (1. - coord_2) +
         this->coord_vtx_[4][1] * (1. - coord_0) *
             (1. - coord_1) * coord_2 +
         this->coord_vtx_[5][1] * coord_0 * (1. - coord_1) *
             coord_2 +
         this->coord_vtx_[6][1] * coord_0 * coord_1 * coord_2 +
         this->coord_vtx_[7][1] * (1. - coord_0) * coord_1 *
             coord_2;
}

template < class DataType, int DIM >
DataType TriLinearHexahedronTransformation< DataType, DIM >::y_x( const Coord &coord_ref) const {
  assert (DIM == 3);
  const DataType coord_1 = coord_ref[1];
  const DataType coord_2 = coord_ref[2];

  return -this->coord_vtx_[0][1] * (1. - coord_1) *
             (1. - coord_2) +
         this->coord_vtx_[1][1] * (1. - coord_1) *
             (1. - coord_2) +
         this->coord_vtx_[2][1] * coord_1 * (1. - coord_2) -
         this->coord_vtx_[3][1] * coord_1 * (1. - coord_2) -
         this->coord_vtx_[4][1] * (1. - coord_1) * coord_2 +
         this->coord_vtx_[5][1] * (1. - coord_1) * coord_2 +
         this->coord_vtx_[6][1] * coord_1 * coord_2 -
         this->coord_vtx_[7][1] * coord_1 * coord_2;
}

template < class DataType, int DIM >
DataType TriLinearHexahedronTransformation< DataType, DIM >::y_y( const Coord &coord_ref) const {
  assert (DIM == 3);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_2 = coord_ref[2];

  return -this->coord_vtx_[0][1] * (1. - coord_0) *
             (1. - coord_2) -
         this->coord_vtx_[1][1] * coord_0 * (1. - coord_2) +
         this->coord_vtx_[2][1] * coord_0 * (1. - coord_2) +
         this->coord_vtx_[3][1] * (1. - coord_0) *
             (1. - coord_2) -
         this->coord_vtx_[4][1] * (1. - coord_0) * coord_2 -
         this->coord_vtx_[5][1] * coord_0 * coord_2 +
         this->coord_vtx_[6][1] * coord_0 * coord_2 +
         this->coord_vtx_[7][1] * (1. - coord_0) * coord_2;
}

template < class DataType, int DIM >
DataType TriLinearHexahedronTransformation< DataType, DIM >::y_z( const Coord &coord_ref) const {
  assert (DIM == 3);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_1 = coord_ref[1];

  return -this->coord_vtx_[0][1] * (1. - coord_0) *
             (1. - coord_1) -
         this->coord_vtx_[1][1] * coord_0 * (1. - coord_1) -
         this->coord_vtx_[2][1] * coord_0 * coord_1 -
         this->coord_vtx_[3][1] * (1. - coord_0) * coord_1 +
         this->coord_vtx_[4][1] * (1. - coord_0) *
             (1. - coord_1) +
         this->coord_vtx_[5][1] * coord_0 * (1. - coord_1) +
         this->coord_vtx_[6][1] * coord_0 * coord_1 +
         this->coord_vtx_[7][1] * (1. - coord_0) * coord_1;
}

template < class DataType, int DIM >
DataType TriLinearHexahedronTransformation< DataType, DIM >::y_xy( const Coord &coord_ref) const {
  assert (DIM == 3);
  const DataType coord_2 = coord_ref[2];

  return +this->coord_vtx_[0][1] * (1. - coord_2) -
         this->coord_vtx_[1][1] * (1. - coord_2) +
         this->coord_vtx_[2][1] * (1. - coord_2) -
         this->coord_vtx_[3][1] * (1. - coord_2) +
         this->coord_vtx_[4][1] * coord_2 -
         this->coord_vtx_[5][1] * coord_2 +
         this->coord_vtx_[6][1] * coord_2 -
         this->coord_vtx_[7][1] * coord_2;
}

template < class DataType, int DIM >
DataType TriLinearHexahedronTransformation< DataType, DIM >::y_xz( const Coord &coord_ref) const {
  assert (DIM == 3);
  const DataType coord_1 = coord_ref[1];

  return +this->coord_vtx_[0][1] * (1. - coord_1) -
         this->coord_vtx_[1][1] * (1. - coord_1) -
         this->coord_vtx_[2][1] * coord_1 +
         this->coord_vtx_[3][1] * coord_1 -
         this->coord_vtx_[4][1] * (1. - coord_1) +
         this->coord_vtx_[5][1] * (1. - coord_1) +
         this->coord_vtx_[6][1] * coord_1 -
         this->coord_vtx_[7][1] * coord_1;
}

template < class DataType, int DIM >
DataType TriLinearHexahedronTransformation< DataType, DIM >::y_yz( const Coord &coord_ref) const {
  assert (DIM == 3);
  const DataType coord_0 = coord_ref[0];

  return +this->coord_vtx_[0][1] * (1. - coord_0) +
         this->coord_vtx_[1][1] * coord_0 -
         this->coord_vtx_[2][1] * coord_0 -
         this->coord_vtx_[3][1] * (1. - coord_0) -
         this->coord_vtx_[4][1] * (1. - coord_0) -
         this->coord_vtx_[5][1] * coord_0 +
         this->coord_vtx_[6][1] * coord_0 +
         this->coord_vtx_[7][1] * (1. - coord_0);
}

template < class DataType, int DIM >
DataType TriLinearHexahedronTransformation< DataType, DIM >::y_xyz( const Coord &coord_ref) const {
  assert (DIM == 3);

  return - this->coord_vtx_[0][1]
         + this->coord_vtx_[1][1] 
         - this->coord_vtx_[2][1] 
         + this->coord_vtx_[3][1] 
         + this->coord_vtx_[4][1] 
         - this->coord_vtx_[5][1] 
         + this->coord_vtx_[6][1]
         - this->coord_vtx_[7][1];
}

template < class DataType, int DIM >
DataType
TriLinearHexahedronTransformation< DataType, DIM >::z(const Coord &coord_ref) const {
  assert (DIM == 3);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_1 = coord_ref[1];
  const DataType coord_2 = coord_ref[2];

  return +this->coord_vtx_[0][2] * (1. - coord_0) *
             (1. - coord_1) * (1. - coord_2) +
         this->coord_vtx_[1][2] * coord_0 * (1. - coord_1) *
             (1. - coord_2) +
         this->coord_vtx_[2][2] * coord_0 * coord_1 *
             (1. - coord_2) +
         this->coord_vtx_[3][2] * (1. - coord_0) * coord_1 *
             (1. - coord_2) +
         this->coord_vtx_[4][2] * (1. - coord_0) *
             (1. - coord_1) * coord_2 +
         this->coord_vtx_[5][2] * coord_0 * (1. - coord_1) *
             coord_2 +
         this->coord_vtx_[6][2] * coord_0 * coord_1 * coord_2 +
         this->coord_vtx_[7][2] * (1. - coord_0) * coord_1 *
             coord_2;
}

template < class DataType, int DIM >
DataType TriLinearHexahedronTransformation< DataType, DIM >::z_x( const Coord &coord_ref) const {
  assert (DIM == 3);
  const DataType coord_1 = coord_ref[1];
  const DataType coord_2 = coord_ref[2];

  return -this->coord_vtx_[0][2] * (1. - coord_1) *
             (1. - coord_2) +
         this->coord_vtx_[1][2] * (1. - coord_1) *
             (1. - coord_2) +
         this->coord_vtx_[2][2] * coord_1 * (1. - coord_2) -
         this->coord_vtx_[3][2] * coord_1 * (1. - coord_2) -
         this->coord_vtx_[4][2] * (1. - coord_1) * coord_2 +
         this->coord_vtx_[5][2] * (1. - coord_1) * coord_2 +
         this->coord_vtx_[6][2] * coord_1 * coord_2 -
         this->coord_vtx_[7][2] * coord_1 * coord_2;
}

template < class DataType, int DIM >
DataType TriLinearHexahedronTransformation< DataType, DIM >::z_y(  const Coord &coord_ref) const {
  assert (DIM == 3);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_2 = coord_ref[2];

  return -this->coord_vtx_[0][2] * (1. - coord_0) *
             (1. - coord_2) -
         this->coord_vtx_[1][2] * coord_0 * (1. - coord_2) +
         this->coord_vtx_[2][2] * coord_0 * (1. - coord_2) +
         this->coord_vtx_[3][2] * (1. - coord_0) *
             (1. - coord_2) -
         this->coord_vtx_[4][2] * (1. - coord_0) * coord_2 -
         this->coord_vtx_[5][2] * coord_0 * coord_2 +
         this->coord_vtx_[6][2] * coord_0 * coord_2 +
         this->coord_vtx_[7][2] * (1. - coord_0) * coord_2;
}

template < class DataType, int DIM >
DataType TriLinearHexahedronTransformation< DataType, DIM >::z_z( const Coord &coord_ref) const {
  assert (DIM == 3);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_1 = coord_ref[1];

  return -this->coord_vtx_[0][2] * (1. - coord_0) *
             (1. - coord_1) -
         this->coord_vtx_[1][2] * coord_0 * (1. - coord_1) -
         this->coord_vtx_[2][2] * coord_0 * coord_1 -
         this->coord_vtx_[3][2] * (1. - coord_0) * coord_1 +
         this->coord_vtx_[4][2] * (1. - coord_0) *
             (1. - coord_1) +
         this->coord_vtx_[5][2] * coord_0 * (1. - coord_1) +
         this->coord_vtx_[6][2] * coord_0 * coord_1 +
         this->coord_vtx_[7][2] * (1. - coord_0) * coord_1;
}

template < class DataType, int DIM >
DataType TriLinearHexahedronTransformation< DataType, DIM >::z_xy( const Coord &coord_ref) const {
  assert (DIM == 3);
  const DataType coord_2 = coord_ref[2];

  return +this->coord_vtx_[0][2] * (1. - coord_2) -
         this->coord_vtx_[1][2] * (1. - coord_2) +
         this->coord_vtx_[2][2] * (1. - coord_2) -
         this->coord_vtx_[3][2] * (1. - coord_2) +
         this->coord_vtx_[4][2] * coord_2 -
         this->coord_vtx_[5][2] * coord_2 +
         this->coord_vtx_[6][2] * coord_2 -
         this->coord_vtx_[7][2] * coord_2;
}

template < class DataType, int DIM >
DataType TriLinearHexahedronTransformation< DataType, DIM >::z_xz( const Coord &coord_ref) const {
  assert (DIM == 3);
  const DataType coord_1 = coord_ref[1];

  return +this->coord_vtx_[0][2] * (1. - coord_1) -
         this->coord_vtx_[1][2] * (1. - coord_1) -
         this->coord_vtx_[2][2] * coord_1 +
         this->coord_vtx_[3][2] * coord_1 -
         this->coord_vtx_[4][2] * (1. - coord_1) +
         this->coord_vtx_[5][2] * (1. - coord_1) +
         this->coord_vtx_[6][2] * coord_1 -
         this->coord_vtx_[7][2] * coord_1;
}

template < class DataType, int DIM >
DataType TriLinearHexahedronTransformation< DataType, DIM >::z_yz( const Coord &coord_ref) const {
  assert (DIM == 3);
  const DataType coord_0 = coord_ref[0];

  return +this->coord_vtx_[0][2] * (1. - coord_0) +
         this->coord_vtx_[1][2] * coord_0 -
         this->coord_vtx_[2][2] * coord_0 -
         this->coord_vtx_[3][2] * (1. - coord_0) -
         this->coord_vtx_[4][2] * (1. - coord_0) -
         this->coord_vtx_[5][2] * coord_0 +
         this->coord_vtx_[6][2] * coord_0 +
         this->coord_vtx_[7][2] * (1. - coord_0);
}

template < class DataType, int DIM >
DataType TriLinearHexahedronTransformation< DataType, DIM >::z_xyz( const Coord &coord_ref) const {
  assert (DIM == 3);

  return - this->coord_vtx_[0][2]
         + this->coord_vtx_[1][2] 
         - this->coord_vtx_[2][2] 
         + this->coord_vtx_[3][2] 
         + this->coord_vtx_[4][2] 
         - this->coord_vtx_[5][2] 
         + this->coord_vtx_[6][2]
         - this->coord_vtx_[7][2];
}

} // namespace doffem
} // namespace hiflow

#endif
