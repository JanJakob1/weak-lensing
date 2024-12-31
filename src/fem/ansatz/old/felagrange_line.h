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

#ifndef __FEM_FELAGRANGE_LINE_H_
#define __FEM_FELAGRANGE_LINE_H_

#include "felagrange.h"
#include <cassert>

namespace hiflow {
namespace doffem {

///
/// \class FELagrangeLine felagrange_line.h
/// \brief Lagrangian Finite Element on a line
/// \author Michael Schick<br>Martin Baumann<br>Julian Kraemer
///

template < class DataType, int DIM >
class FELagrangeLine : public FELagrange< DataType, DIM > {
public:
  typedef Vec<DIM, DataType> Coord;

  /// Default Constructor
  FELagrangeLine();

  /// Default Destructor
  ~FELagrangeLine();

  std::string get_name() const { return "LagrangeLine"; }

  void N(const Coord &pt, std::vector< DataType > &weight) const;

  void N_x(const Coord &pt, std::vector< DataType > &weight) const;
  void N_y(const Coord &pt, std::vector< DataType > &weight) const;
  void N_z(const Coord &pt, std::vector< DataType > &weight) const;

  void N_xx(const Coord &pt, std::vector< DataType > &weight) const;
  void N_xy(const Coord &pt, std::vector< DataType > &weight) const;
  void N_xz(const Coord &pt, std::vector< DataType > &weight) const;
  void N_yy(const Coord &pt, std::vector< DataType > &weight) const;
  void N_yz(const Coord &pt, std::vector< DataType > &weight) const;
  void N_zz(const Coord &pt, std::vector< DataType > &weight) const;

protected:
  void init_coord();

private:
  /// Index ordering from plane (x=i, y=j) to vector index
  int ij2ind(int i, int j) const;
};

template < class DataType, int DIM > FELagrangeLine< DataType, DIM >::FELagrangeLine() {
  this->my_type_ = FEType< DataType, DIM >::LAGRANGE_LINE;

  // initialize reference cell

  assert(this->ref_cell_ == NULL);
  this->ref_cell_ = &(mesh::CellType::get_instance(mesh::CellType::LINE));
}

template < class DataType, int DIM > FELagrangeLine< DataType, DIM >::~FELagrangeLine() {}

template < class DataType, int DIM > void FELagrangeLine< DataType, DIM >::init_coord() {
  assert (DIM == 1);
  // set topological degree

  this->tdim_ = 1;

  if (this->my_deg_ == 0) {
    this->coord_.clear();

    // Coordinates of the middle-point of line
    Coord coord;

    // Centroid
    coord[0] = 0.5;

    this->coord_.push_back(coord);
  } else {
    assert(this->my_deg_ > 0);

    // Lexicographical ordering

    const DataType offset = (1.0 / this->my_deg_);

    this->coord_.clear();
    const int nb_dof_on_cell = (this->my_deg_ + 1);
    this->coord_.resize(nb_dof_on_cell);

    const int nb_dof_line = this->my_deg_ + 1;

    // Filling coord vector by lexicographical strategy

    for (int j = 0; j < nb_dof_line; ++j) {
      Coord coord;

      coord[0] = j * offset;

      this->coord_[ij2ind(0, j)] = coord;
    }
  }
}

template < class DataType, int DIM >
int FELagrangeLine< DataType, DIM >::ij2ind(int i, int j) const {
  assert (DIM == 1);

  return (j);
}

/// \details Every degree of a used lagrangian polynomial has to satisfy the
///          condition degree < this->my_deg_. All possible combinations are
///          multiplied and considered in the sum.

template < class DataType, int DIM >
void FELagrangeLine< DataType, DIM >::N(const Coord &pt,
                                   std::vector< DataType > &weight) const {
  assert (DIM == 1);
  assert(weight.size() == this->get_nb_dof_on_cell());

  if (this->my_deg_ > 0) {
    for (int j = 0; j <= this->my_deg_; ++j) {
      weight[ij2ind(0, j)] = this->lp_.poly(this->my_deg_, j, pt[0]);
    }
  } else {
    weight[ij2ind(0, 0)] = 1.0;
  }
}

/// \details Every degree of a used lagrangian polynomial has to satisfy the
///          condition degree < this->my_deg_. All possible combinations are
///          multiplied and considered in the sum, w.r.t. the derivatives for
///          the x - variable.

template < class DataType, int DIM >
void FELagrangeLine< DataType, DIM >::N_x(const Coord &pt,
                                     std::vector< DataType > &weight) const {
  assert (DIM == 1);
  assert(weight.size() == this->get_nb_dof_on_cell());

  if (this->my_deg_ > 0) {
    for (int j = 0; j <= this->my_deg_; ++j)

    {
      weight[ij2ind(0, j)] = this->lp_.poly_x(this->my_deg_, j, pt[0]);
    }
  } else {
    weight[ij2ind(0, 0)] = 0.0;
  }
}

/// \details There are no y - derivatives in a 1D case

template < class DataType, int DIM >
void FELagrangeLine< DataType, DIM >::N_y(const Coord &pt,
                                     std::vector< DataType > &weight) const {
  assert (DIM == 1);
  weight.assign(weight.size(), 0.);
}

/// \details There are no z - derivatives in a 1D case

template < class DataType, int DIM >
void FELagrangeLine< DataType, DIM >::N_z(const Coord &pt,
                                     std::vector< DataType > &weight) const {
  assert (DIM == 1);
  weight.assign(weight.size(), 0.);
}

/// \details Every degree of a used lagrangian polynomial has to satisfy the
///          condition degree < this->my_deg_. All possible combinations are
///          multiplied and considered in the sum, w.r.t. the second derivatives
///          for the xx - variable.

template < class DataType, int DIM >
void FELagrangeLine< DataType, DIM >::N_xx(const Coord &pt,
                                      std::vector< DataType > &weight) const {
  assert (DIM == 1);
  assert(weight.size() == this->get_nb_dof_on_cell());

  if (this->my_deg_ > 0) {
    for (int j = 0; j <= this->my_deg_; ++j) {
      weight[ij2ind(0, j)] = this->lp_.poly_xx(this->my_deg_, j, pt[0]);
    }
  } else {
    weight[ij2ind(0, 0)] = 0.0;
  }
}

/// \details There are no y - derivatives in a 1D case

template < class DataType, int DIM >
void FELagrangeLine< DataType, DIM >::N_xy(const Coord &pt,
                                      std::vector< DataType > &weight) const {
  assert (DIM == 1);
  weight.assign(weight.size(), 0.);
}

/// \details There are no z - derivatives in a 1D case

template < class DataType, int DIM >
void FELagrangeLine< DataType, DIM >::N_xz(const Coord &pt,
                                      std::vector< DataType > &weight) const {
  assert (DIM == 1);
  weight.assign(weight.size(), 0.);
}

/// \details There are no y - derivatives in a 1D case

template < class DataType, int DIM >
void FELagrangeLine< DataType, DIM >::N_yy(const Coord &pt,
                                      std::vector< DataType > &weight) const {
  assert (DIM == 1);
  weight.assign(weight.size(), 0.);
}

/// \details There are no z - derivatives in a 1D case

template < class DataType, int DIM >
void FELagrangeLine< DataType, DIM >::N_yz(const Coord &pt,
                                      std::vector< DataType > &weight) const {
  assert (DIM == 1);
  weight.assign(weight.size(), 0.);
}

/// \details There are no z - derivatives in a 1D case

template < class DataType, int DIM >
void FELagrangeLine< DataType, DIM >::N_zz(const Coord &pt,
                                      std::vector< DataType > &weight) const {
  assert (DIM == 1);
  weight.assign(weight.size(), 0.);
}

} // namespace doffem
} // namespace hiflow
#endif
