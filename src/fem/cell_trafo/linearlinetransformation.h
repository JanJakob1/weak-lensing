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

#ifndef __FEM_LINEAR_LINE_TRANSFORMATION_H_
#define __FEM_LINEAR_LINE_TRANSFORMATION_H_

#include "fem/cell_trafo/cell_transformation.h"
#include <cmath>
#include <iomanip>
#include <iostream>

namespace hiflow {
namespace doffem {

///
/// \class LinearLineTransformation linearlinetransformation.h
/// \brief Linear transformation mapping from reference to physical cell for a
/// Line \author Michael Schick<br>Martin Baumann<br>Julian Kraemer
///

template < class DataType, int DIM >
class LinearLineTransformation : public CellTransformation< DataType, DIM > {
public:
  typedef Vec<DIM, DataType> Coord;

  explicit LinearLineTransformation(ConstRefCellPtr<DataType, DIM> ref_cell);

  bool differs_by_translation_from (ConstCellTrafoPtr<DataType, DIM> rhs) const;
   
  bool inverse(Coord co_phy, Coord &co_ref) const;
  bool inverse_1Dto2D(Vec<2, DataType> co_phy, Coord &co_ref) const;
  bool inverse_1Dto3D(Vec<3, DataType> co_phy, Coord &co_ref) const;

  DataType x(const Coord &coord_ref) const;
  DataType x_x(const Coord &coord_ref) const;
  DataType y(const Coord & coord_ref) const;
  DataType y_x(const Coord & coord_ref) const;
  DataType z(const Coord & coord_ref) const;
  DataType z_x(const Coord & coord_ref) const;

};

template < class DataType, int DIM >
LinearLineTransformation< DataType, DIM >::LinearLineTransformation(ConstRefCellPtr<DataType, DIM> ref_cell)
    : CellTransformation< DataType, DIM >(ref_cell) 
{
  this->order_ = 1;
  this->fixed_ref_cell_type_ = REF_CELL_LINE_STD;
  this->name_ = "Line";
}

template < class DataType, int DIM >
bool LinearLineTransformation< DataType, DIM >::differs_by_translation_from (ConstCellTrafoPtr<DataType, DIM> rhs) const
{
  assert(this->ref_cell_);
  
  if (this->ref_cell_->type() != rhs->get_ref_cell()->type())
  {
    return false;
  }
  
  std::vector< Coord > ref_coords = this->ref_cell_->get_coords();
  assert (ref_coords.size() == 2);
  
  Coord my_p01 = this->transform(ref_coords[0]) - this->transform(ref_coords[1]); 
  Coord rhs_p01 = rhs->transform(ref_coords[0]) - rhs->transform(ref_coords[1]); 
  if (my_p01 != rhs_p01)
  {
    return false;
  }
  return true;
}

template < class DataType, int DIM >
bool LinearLineTransformation< DataType, DIM >::inverse(Coord co_phy, Coord &co_ref) const {

  co_ref[0] = (co_phy[0] - this->coord_vtx_[0][0]) / (this->coord_vtx_[1][0] - this->coord_vtx_[0][0]);
  return this->contains_reference_point(co_ref);
}

template < class DataType, int DIM >
DataType LinearLineTransformation< DataType, DIM >::x(const Coord &coord_ref) const {
  return coord_ref[0] * (this->coord_vtx_[1][0] - this->coord_vtx_[0][0]) 
        + this->coord_vtx_[0][0];
}


template < class DataType, int DIM >
DataType
LinearLineTransformation< DataType, DIM >::x_x(const Coord &coord_ref) const {
  return this->coord_vtx_[1][0] - this->coord_vtx_[0][0];
}

// begin preliminary test functions when the dimension of the physical point is > 1

template < class DataType, int DIM >
DataType
LinearLineTransformation< DataType, DIM >:: y(const Coord & coord_ref) const {
  return coord_ref[0] * (this->coord_vtx_[1][1] - this->coord_vtx_[0][1]) 
        + this->coord_vtx_[0][1];
}

template < class DataType, int DIM >
DataType
LinearLineTransformation< DataType, DIM >:: z(const Coord & coord_ref) const {
  return coord_ref[0] * (this->coord_vtx_[1][2] - this->coord_vtx_[0][2]) 
        + this->coord_vtx_[0][2];
}

template < class DataType, int DIM >
DataType
LinearLineTransformation< DataType, DIM >::y_x(const Coord &coord_ref) const {
  return this->coord_vtx_[1][1] - this->coord_vtx_[0][1];
}

template < class DataType, int DIM >
DataType
LinearLineTransformation< DataType, DIM >::z_x(const Coord &coord_ref) const {
  return this->coord_vtx_[1][2] - this->coord_vtx_[0][2];
}

//The inverse function can actually remain the same. Because all the rows in the linear system are pairwise linearly dependent, it suffices to solve one
//However, one can check if the physial point actually lies in the span of the physical cell by checking if the solutions provided by each row individually coincide.

template < class DataType, int DIM >
bool LinearLineTransformation< DataType, DIM >::inverse_1Dto2D(Vec<2, DataType> co_phy, Coord &co_ref) const {
  DataType sol_0 = (co_phy[0] - this->coord_vtx_[0][0]) 
                 / (this->coord_vtx_[1][0] - this->coord_vtx_[0][0]);
  DataType sol_1 = (co_phy[1] - this->coord_vtx_[0][1]) 
                 / (this->coord_vtx_[1][1] - this->coord_vtx_[0][1]);
  
  if (sol_0 != sol_1) 
  {
    return false;
  }
  
  co_ref[0] = sol_0;
  return this->contains_reference_point(co_ref);
}
template < class DataType, int DIM >
bool LinearLineTransformation< DataType, DIM >::inverse_1Dto3D(Vec<3, DataType> co_phy, Coord &co_ref) const {
  DataType sol_0 = (co_phy[0] - this->coord_vtx_[0][0]) 
                 / (this->coord_vtx_[1][0] - this->coord_vtx_[0][0]);
  DataType sol_1 = (co_phy[1] - this->coord_vtx_[0][1]) 
                 / (this->coord_vtx_[1][1] - this->coord_vtx_[0][1]);
  DataType sol_2 = (co_phy[2] - this->coord_vtx_[0][2])
                 / (this->coord_vtx_[1][2] - this->coord_vtx_[0][2]);
  
  if ((sol_0 != sol_1) || (sol_0!= sol_1) || (sol_1 != sol_2) )
  { 
    return false;
  }
  
  co_ref[0] = sol_0;
  return this->contains_reference_point(co_ref);
}


} // namespace doffem
} // namespace hiflow

#endif
