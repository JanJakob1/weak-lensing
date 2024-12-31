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

#ifndef __FEM_BILINEAR_QUAD_TRANSFORMATION_H_
#define __FEM_BILINEAR_QUAD_TRANSFORMATION_H_

#include "fem/cell_trafo/cell_transformation.h"
#include "common/log.h"
#include "fem/cell_trafo/lineartriangletransformation.h"
#include <cmath>
#include <iomanip>

namespace hiflow {
namespace doffem {

///
/// \class BiLinearQuadTransformation bilinearquadtransformation.h
/// \brief Bilinear transformation mapping from reference to physical cell for a
/// Quadrilateral /// \author Michael Schick<br>Martin Baumann<br>Simon Gawlok<br>Philipp Gerstner
///

template < class DataType, int DIM >
class BiLinearQuadTransformation : public CellTransformation< DataType, DIM > {
public:
  typedef Vec<DIM, DataType> Coord;

  explicit BiLinearQuadTransformation(ConstRefCellPtr<DataType, DIM> ref_cell);

  ~ BiLinearQuadTransformation ()
  {
  }

  bool differs_by_translation_from (ConstCellTrafoPtr<DataType, DIM> rhs) const;
   
  bool inverse(Coord co_phy, Coord &co_ref) const;
  
  bool inverse_2Dto3D(Vec<3, DataType> co_phy, Vec<2, DataType> &co_ref) const;

  DataType x(const Coord &coord_ref) const;
  DataType x_x(const Coord &coord_ref) const;
  DataType x_y(const Coord &coord_ref) const;
  DataType x_xy(const Coord &coord_ref) const;

  DataType y(const Coord &coord_ref) const;
  DataType y_x(const Coord &coord_ref) const;
  DataType y_y(const Coord &coord_ref) const;
  DataType y_xy(const Coord &coord_ref) const;
  
  DataType z(const Coord &coord_ref) const;
  DataType z_x(const Coord &coord_ref) const;
  DataType z_y(const Coord &coord_ref) const;
  DataType z_xy(const Coord &coord_ref) const;

protected:
  bool inverse_by_decomposition(Coord co_phy, Coord &co_ref) const;
  
  bool inverse_by_decomposition_2Dto3D(Vec<3, DataType> co_phy, Vec<2, DataType> &co_ref) const; 

  int8_t triangle_decomp_ind_[2][3];
};

template < class DataType, int DIM >
BiLinearQuadTransformation< DataType, DIM >::BiLinearQuadTransformation(ConstRefCellPtr<DataType, DIM> ref_cell)
    : CellTransformation< DataType, DIM >(ref_cell) 
{
  this->order_ = 2;
  this->fixed_ref_cell_type_ = REF_CELL_QUAD_STD;
  this->name_ = "Quad";
    
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

template < class DataType, int DIM >
bool BiLinearQuadTransformation< DataType, DIM >::differs_by_translation_from (ConstCellTrafoPtr<DataType, DIM> rhs) const
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
bool BiLinearQuadTransformation< DataType, DIM >::inverse_by_decomposition(Coord co_phy, Coord &co_ref) const 
{
  assert (DIM == 2);
  RefCellPtr <DataType, DIM> ref_tri = RefCellPtr<DataType, DIM>(new RefCellTriStd<DataType, DIM>);
  LinearTriangleTransformation< DataType, DIM > tria_trafo(ref_tri);
  Coord co_ref_tri;

  int tria_ind = -1;

  // loop through tetrahedron decomposition
  for (int t = 0; t < 2; ++t) 
  {
    // build linear tetrahedron transformation
    std::vector<DataType> tria_coord(6);
    for (int v = 0; v < 3; ++v) 
    {
      for (int d = 0; d < 2; ++d) 
      {
        tria_coord[this->ij2ind(v,d)] = this->coord_vtx_[this->triangle_decomp_ind_[t][v]][d];
      }
    }
    tria_trafo.reinit(tria_coord);

    // compute reference coordinates w.r.t. tetrahedron
    bool tria_success = tria_trafo.inverse(co_phy, co_ref_tri);
    if (!tria_success) 
    {
      continue;
    }

    // check whether reference coordinates are contained in tetrahedron

    // if yes: convert reference coordinates
    if (tria_trafo.contains_reference_point(co_ref_tri)) 
    {
      switch (t) 
      {
      case 0:
        co_ref = co_ref_tri;
        break;
      case 1:
        co_ref[0] = 1. - co_ref_tri[0];
        co_ref[1] = 1. - co_ref_tri[1];
        break;
      }
      tria_ind = t;
      break;
    }
  }
  // TODO: is this check necessary?
  // check if reference point found by decomposition is sufficiently accurate
  Coord phy_test;
  phy_test[0] = this->x(co_ref);
  phy_test[1] = this->y(co_ref);

  DataType residual = norm(phy_test - co_phy);

#ifndef NDEBUG
  LOG_DEBUG(2, "==============");
  LOG_DEBUG(3, "Physical point " << co_phy[0] << " " << co_phy[1]
                                 << " is contained in triangle " << tria_ind);

  if (tria_ind == -1) {
    LOG_DEBUG(2, "Physical point " << co_phy[0] << " " << co_phy[1]
                                   << " is contained in triangle " << tria_ind);
    for (int v = 0; v < 4; ++v) {
      LOG_DEBUG(2, "vertex " << v << " : "
                             << this->coord_vtx_[v][0] << ", "
                             << this->coord_vtx_[v][1]);
    }
  }

  LOG_DEBUG(
      3, "Inverse cell-trafo based on triangle decomposition yields residual = "
             << residual);
#endif

  if (residual < 5. * std::numeric_limits< DataType >::epsilon()) 
  {
    return true;
  } 
  else 
  {
    return false;
  }
}

template < class DataType, int DIM >
bool BiLinearQuadTransformation< DataType, DIM >::inverse(Coord co_phy, Coord &co_ref) const {
  assert (DIM == 2);
  bool success = this->inverse_by_decomposition(co_phy, co_ref);

  if (success) 
  {
    // reference point obtained by decomposition is sufficiently accurate ->
    // exit
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
    /*
    // Newton method with old initial value
    if (!success)
    {
        success = this->inverse_newton_2d ( x_phy, y_phy, x_ref, y_ref );
    }
    * */
  }
}

template < class DataType, int DIM >
bool BiLinearQuadTransformation< DataType, DIM >::inverse_by_decomposition_2Dto3D(Vec<3, DataType> co_phy, Vec<2, DataType> &co_ref) const 
{
  ConstRefCellPtr< DataType, DIM > ref_tri = ConstRefCellPtr< DataType, DIM >(new RefCellTriStd<DataType, DIM>);
  LinearTriangleTransformation< DataType, DIM > tria_trafo(ref_tri);
  Vec<2, DataType> co_ref_tri;


  int tria_ind = -1;

  // loop through tetrahedron decomposition
  for (int t = 0; t < 2; ++t) 
  {
    // build linear tetrahedron transformation
    std::vector<DataType> tria_coord(9);
    for (int v = 0; v < 3; ++v) 
    {
      for (int d = 0; d < 3; ++d) 
      {
        tria_coord[this->ij2ind(v,d)] = this->coord_vtx_[this->triangle_decomp_ind_[t][v]][d];    //TODO: coord_vtx_ anpassen
      }
    }
    tria_trafo.reinit(tria_coord);

    // compute reference coordinates w.r.t. tetrahedron
    // TODO solve template mismatch problem
    //bool tria_success = tria_trafo.inverse(co_phy, co_ref_tri);
    bool tria_success = false;
    if (!tria_success) 
    {
      continue;
    }

    // check whether reference coordinates are contained in tetrahedron

    // if yes: convert reference coordinates
    if (tria_trafo.contains_reference_point(co_ref_tri)) 
    {
      switch (t) 
      {
      case 0:
        co_ref = co_ref_tri;
        break;
      case 1:
        co_ref[0] = 1. - co_ref_tri[0];
        co_ref[1] = 1. - co_ref_tri[1];
        break;
      }
      tria_ind = t;
      break;
    }
  }
  // TODO: is this check necessary?
  // check if reference point found by decomposition is sufficiently accurate
  Vec<3, DataType> phy_test;
  phy_test[0] = this->x(co_ref);
  phy_test[1] = this->y(co_ref);
  phy_test[2] = this->z(co_ref);
 
  DataType residual = norm(phy_test - co_phy);

#ifndef NDEBUG
  LOG_DEBUG(2, "==============");
  LOG_DEBUG(3, "Physical point " << co_phy[0] << " " << co_phy[1]
                                 << " is contained in triangle " << tria_ind);

  if (tria_ind == -1) {
    LOG_DEBUG(2, "Physical point " << co_phy[0] << " " << co_phy[1]
                                   << " is contained in triangle " << tria_ind);
    for (int v = 0; v < 4; ++v) {
      LOG_DEBUG(2, "vertex " << v << " : "
                             << this->coord_vtx_[v][0] << ", "
                             << this->coord_vtx_[v][1]);
    }
  }

  LOG_DEBUG(
      3, "Inverse cell-trafo based on triangle decomposition yields residual = "
             << residual);
#endif

  if (residual < 5. * std::numeric_limits< DataType >::epsilon()) 
  {
    return true;
  } 
  else 
  {
    return false;
  }
}

template < class DataType, int DIM >
bool BiLinearQuadTransformation< DataType, DIM >::inverse_2Dto3D(Vec<3, DataType> co_phy, Vec<2, DataType> &co_ref) const {
 
  bool success = this->inverse_by_decomposition_2Dto3D(co_phy, co_ref);

  if (success) 
  {
    // reference point obtained by decomposition is sufficiently accurate ->
    // exit
    return true;
  } 
  else 
  {
    // reference point obtained by decomposition is not accurate enough -> use
    // as initial value for Newton method
    bool newton_success = this->inverse_newton_2Dto3D(co_phy, co_ref, co_ref);
    if (newton_success)
    {
      return this->contains_reference_point(co_ref);
    }
    else
    {
      return false;
    }
    /*
    // Newton method with old initial value
    if (!success)
    {
        success = this->inverse_newton_2d ( x_phy, y_phy, x_ref, y_ref );
    }
    * */
  }
}

template < class DataType, int DIM >
DataType
BiLinearQuadTransformation< DataType, DIM >::x(const Coord &coord_ref) const {
  assert (DIM == 2);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_1 = coord_ref[1];

  return this->coord_vtx_[0][0] * (1. - coord_0) * (1. - coord_1) +
         this->coord_vtx_[1][0] * coord_0 * (1. - coord_1) +
         this->coord_vtx_[2][0] * coord_0 * coord_1 +
         this->coord_vtx_[3][0] * (1. - coord_0) * coord_1;
}

template < class DataType, int DIM >
DataType BiLinearQuadTransformation< DataType, DIM >::x_x(const Coord &coord_ref) const {
  assert (DIM == 2);
  const DataType coord_1 = coord_ref[1];

  return this->coord_vtx_[0][0] * (-1. + coord_1) +
         this->coord_vtx_[1][0] * (1. - coord_1) +
         this->coord_vtx_[2][0] * coord_1 -
         this->coord_vtx_[3][0] * coord_1;
}

template < class DataType, int DIM >
DataType BiLinearQuadTransformation< DataType, DIM >::x_y(const Coord &coord_ref) const {
  assert (DIM == 2);
  const DataType coord_0 = coord_ref[0];

  return this->coord_vtx_[0][0] * (-1. + coord_0) -
         this->coord_vtx_[1][0] * coord_0 +
         this->coord_vtx_[2][0] * coord_0 +
         this->coord_vtx_[3][0] * (1. - coord_0);
}

template < class DataType, int DIM >
DataType BiLinearQuadTransformation< DataType, DIM >::x_xy(const Coord &coord_ref) const {
  assert (DIM == 2);
  return this->coord_vtx_[0][0] -
         this->coord_vtx_[1][0] +
         this->coord_vtx_[2][0] -
         this->coord_vtx_[3][0];
}

template < class DataType, int DIM >
DataType BiLinearQuadTransformation< DataType, DIM >::y(const Coord &coord_ref) const {
  assert (DIM == 2);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_1 = coord_ref[1];

  return this->coord_vtx_[0][1] * (1. - coord_0) * (1. - coord_1) +
         this->coord_vtx_[1][1] * coord_0 * (1. - coord_1) +
         this->coord_vtx_[2][1] * coord_0 * coord_1 +
         this->coord_vtx_[3][1] * (1. - coord_0) * coord_1;
}

template < class DataType, int DIM >
DataType BiLinearQuadTransformation< DataType, DIM >::y_x(const Coord &coord_ref) const {
  assert (DIM == 2);
  const DataType coord_1 = coord_ref[1];

  return this->coord_vtx_[0][1] * (-1. + coord_1) +
         this->coord_vtx_[1][1] * (1. - coord_1) +
         this->coord_vtx_[2][1] * coord_1 -
         this->coord_vtx_[3][1] * coord_1;
}

template < class DataType, int DIM >
DataType
BiLinearQuadTransformation< DataType, DIM >::y_y(const Coord &coord_ref) const {
  assert (DIM == 2);
  const DataType coord_0 = coord_ref[0];

  return this->coord_vtx_[0][1] * (-1. + coord_0) -
         this->coord_vtx_[1][1] * coord_0 +
         this->coord_vtx_[2][1] * coord_0 +
         this->coord_vtx_[3][1] * (1. - coord_0);
}

template < class DataType, int DIM >
DataType BiLinearQuadTransformation< DataType, DIM >::y_xy(const Coord &coord_ref) const {
  assert (DIM == 2);
  return this->coord_vtx_[0][1] -
         this->coord_vtx_[1][1] +
         this->coord_vtx_[2][1] -
         this->coord_vtx_[3][1];
}

// begin preliminary test functions when the dimension of the physical point is > 1
template < class DataType, int DIM >
DataType BiLinearQuadTransformation< DataType, DIM >::z(const Coord &coord_ref) const {
  assert (DIM == 2);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_1 = coord_ref[1];

  return this->coord_vtx_[0][2] * (1. - coord_0) * (1. - coord_1) +
         this->coord_vtx_[1][2] * coord_0 * (1. - coord_1) +
         this->coord_vtx_[2][2] * coord_0 * coord_1 +
         this->coord_vtx_[3][2] * (1. - coord_0) * coord_1;
}

template < class DataType, int DIM >
DataType BiLinearQuadTransformation< DataType, DIM >::z_x(const Coord &coord_ref) const {
  assert (DIM == 2);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_1 = coord_ref[1];

  return this->coord_vtx_[0][2] * (coord_1 - 1.) +
         this->coord_vtx_[1][2] * (1. - coord_1) +
         this->coord_vtx_[2][2] * coord_1 -
         this->coord_vtx_[3][2] * coord_1;
}

template < class DataType, int DIM >
DataType BiLinearQuadTransformation< DataType, DIM >::z_y(const Coord &coord_ref) const {
  assert (DIM == 2);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_1 = coord_ref[1];

  return this->coord_vtx_[0][2] * (coord_0 - 1.) -
         this->coord_vtx_[1][2] * coord_0  +
         this->coord_vtx_[2][2] * coord_0 +
         this->coord_vtx_[3][2] * (1. - coord_0);
}

template < class DataType, int DIM >
DataType BiLinearQuadTransformation< DataType, DIM >::z_xy(const Coord &coord_ref) const {
  assert (DIM == 2);
  const DataType coord_0 = coord_ref[0];
  const DataType coord_1 = coord_ref[1];

  return this->coord_vtx_[0][2] -
         this->coord_vtx_[1][2] +
         this->coord_vtx_[2][2] -
         this->coord_vtx_[3][2];
}


} // namespace doffem
} // namespace hiflow

#endif
