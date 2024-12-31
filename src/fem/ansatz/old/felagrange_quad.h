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

#ifndef __FEM_FELAGRANGE_QUAD_H_
#define __FEM_FELAGRANGE_QUAD_H_

#include "felagrange.h"
#include <cassert>

namespace hiflow {
namespace doffem {

///
/// \class FELagrangeQuad felagrange_quad.h
/// \brief Lagrangian Finite Element on a Quadrilateral
/// \author Michael Schick<br>Martin Baumann
///

template < class DataType, int DIM >
class FELagrangeQuad : public FELagrange< DataType, DIM > {
public:
  typedef Vec<DIM, DataType> Coord;

  /// Default Constructor
  FELagrangeQuad();

  /// Default Destructor
  ~FELagrangeQuad();

  std::string get_name() const { return "LagrangeQuadrilateral"; }

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

template < class DataType, int DIM > FELagrangeQuad< DataType, DIM >::FELagrangeQuad() {
  this->my_type_ = FEType< DataType, DIM >::LAGRANGE_QUAD;

  // initialize reference cell

  assert(this->ref_cell_ == NULL);
  this->ref_cell_ =
      &(mesh::CellType::get_instance(mesh::CellType::QUADRILATERAL));
}

template < class DataType, int DIM > FELagrangeQuad< DataType, DIM >::~FELagrangeQuad() {}

template < class DataType, int DIM > void FELagrangeQuad< DataType, DIM >::init_coord() {
  // set topological degree
  assert (DIM == 2);
  this->tdim_ = 2;

  if (this->my_deg_ == 0) {
    this->coord_.clear();

    // Coordinates of the middle-point of quad
    Coord coord;

    // Centroid
    coord[0] = 0.5;
    coord[1] = 0.5;

    this->coord_.push_back(coord);
  } else {
    assert(this->my_deg_ > 0);

    // Lexicographical ordering

    const DataType offset = (1.0 / this->my_deg_);

    this->coord_.clear();
    const int nb_dof_on_cell = (this->my_deg_ + 1) * (this->my_deg_ + 1);
    this->coord_.resize(nb_dof_on_cell);

    const int nb_dof_line = this->my_deg_ + 1;

    assert(this->coord_.size() == nb_dof_line * nb_dof_line);

    // Filling coord vector by lexicographical strategy

    for (int j = 0; j < nb_dof_line; ++j) {
      const DataType j_offset = j * offset;
      for (int i = 0; i < nb_dof_line; ++i) {
        Coord coord;

        coord[0] = i * offset;
        coord[1] = j_offset;

        this->coord_[ij2ind(i, j)] = coord;
      }
    }
  }
}

template < class DataType, int DIM >
int FELagrangeQuad< DataType, DIM >::ij2ind(int i, int j) const {
  assert (DIM == 2);
  const int nb_dof_on_line = this->my_deg_ + 1;
  return (i + j * nb_dof_on_line);
}

/// \details Every degree of a used lagrangian polynomial has to satisfy the
///          condition degree < this->my_deg_. All possible combinations are
///          multiplied and considered in the sum.

template < class DataType, int DIM >
void FELagrangeQuad< DataType, DIM >::N(const Coord &pt,
                                   std::vector< DataType > &weight) const {
  assert (DIM == 2);
  assert(weight.size() == this->get_nb_dof_on_cell());

  if (this->my_deg_ > 0) {
    for (int j = 0; j <= this->my_deg_; ++j) {
      const DataType lp_j = this->lp_.poly(this->my_deg_, j, pt[1]);
      for (int i = 0; i <= this->my_deg_; ++i) {
        weight[ij2ind(i, j)] = this->lp_.poly(this->my_deg_, i, pt[0]) * lp_j;
      }
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
void FELagrangeQuad< DataType, DIM >::N_x(const Coord &pt,
                                     std::vector< DataType > &weight) const {
  assert (DIM == 2);
  assert(weight.size() == this->get_nb_dof_on_cell());

  if (this->my_deg_ > 0) {
    for (int j = 0; j <= this->my_deg_; ++j) {
      const DataType lp_j = this->lp_.poly(this->my_deg_, j, pt[1]);
      for (int i = 0; i <= this->my_deg_; ++i) {
        weight[ij2ind(i, j)] = this->lp_.poly_x(this->my_deg_, i, pt[0]) * lp_j;
      }
    }
  } else {
    weight[ij2ind(0, 0)] = 0.0;
  }
}

/// \details Every degree of a used lagrangian polynomial has to satisfy the
///          condition degree < this->my_deg_. All possible combinations are
///          multiplied and considered in the sum, w.r.t. the derivatives for
///          the y - variable.

template < class DataType, int DIM >
void FELagrangeQuad< DataType, DIM >::N_y(const Coord &pt,
                                     std::vector< DataType > &weight) const {
  assert (DIM == 2);
  assert(weight.size() == this->get_nb_dof_on_cell());

  if (this->my_deg_ > 0) {
    for (int j = 0; j <= this->my_deg_; ++j) {
      const DataType lp_j = this->lp_.poly_x(this->my_deg_, j, pt[1]);
      for (int i = 0; i <= this->my_deg_; ++i) {
        weight[ij2ind(i, j)] = this->lp_.poly(this->my_deg_, i, pt[0]) * lp_j;
      }
    }
  } else {
    weight[ij2ind(0, 0)] = 0.0;
  }
}

/// \details There are no z - derivatives in a 2D case

template < class DataType, int DIM >
void FELagrangeQuad< DataType, DIM >::N_z(const Coord &pt,
                                     std::vector< DataType > &weight) const {
  assert (DIM == 2);
  weight.assign(weight.size(), 0.);
}

/// \details Every degree of a used lagrangian polynomial has to satisfy the
///          condition degree < this->my_deg_. All possible combinations are
///          multiplied and considered in the sum, w.r.t. the second derivatives
///          for the xx - variable.

template < class DataType, int DIM >
void FELagrangeQuad< DataType, DIM >::N_xx(const Coord &pt,
                                      std::vector< DataType > &weight) const {
  assert (DIM == 2);
  assert(weight.size() == this->get_nb_dof_on_cell());

  if (this->my_deg_ > 0) {
    for (int j = 0; j <= this->my_deg_; ++j) {
      const DataType lp_j = this->lp_.poly(this->my_deg_, j, pt[1]);
      for (int i = 0; i <= this->my_deg_; ++i) {
        weight[ij2ind(i, j)] =
            this->lp_.poly_xx(this->my_deg_, i, pt[0]) * lp_j;
      }
    }
  } else {
    weight[ij2ind(0, 0)] = 0.0;
  }
}

/// \details Every degree of a used lagrangian polynomial has to satisfy the
///          condition degree < this->my_deg_. All possible combinations are
///          multiplied and considered in the sum, w.r.t. the second derivatives
///          for the xy - variable.

template < class DataType, int DIM >
void FELagrangeQuad< DataType, DIM >::N_xy(const Coord &pt,
                                      std::vector< DataType > &weight) const {
  assert (DIM == 2);
  assert(weight.size() == this->get_nb_dof_on_cell());

  if (this->my_deg_ > 0) {
    for (int j = 0; j <= this->my_deg_; ++j) {
      const DataType lp_j = this->lp_.poly_x(this->my_deg_, j, pt[1]);
      for (int i = 0; i <= this->my_deg_; ++i) {
        weight[ij2ind(i, j)] = this->lp_.poly_x(this->my_deg_, i, pt[0]) * lp_j;
      }
    }
  } else {
    weight[ij2ind(0, 0)] = 0.0;
  }
}

/// \details There are no z - derivatives in a 2D case

template < class DataType, int DIM >
void FELagrangeQuad< DataType, DIM >::N_xz(const Coord &pt,
                                      std::vector< DataType > &weight) const {
  assert (DIM == 2);
  weight.assign(weight.size(), 0.);
}

/// \details Every degree of a used lagrangian polynomial has to satisfy the
///          condition degree < this->my_deg_. All possible combinations are
///          multiplied and considered in the sum, w.r.t. the second derivatives
///          for the yy - variable.

template < class DataType, int DIM >
void FELagrangeQuad< DataType, DIM >::N_yy(const Coord &pt,
                                      std::vector< DataType > &weight) const {
  assert (DIM == 2);
  assert(weight.size() == this->get_nb_dof_on_cell());

  if (this->my_deg_ > 0) {
    for (int j = 0; j <= this->my_deg_; ++j) {
      const double lp_j = this->lp_.poly_xx(this->my_deg_, j, pt[1]);
      for (int i = 0; i <= this->my_deg_; ++i) {
        weight[ij2ind(i, j)] = this->lp_.poly(this->my_deg_, i, pt[0]) * lp_j;
      }
    }
  } else {
    weight[ij2ind(0, 0)] = 0.0;
  }
}

/// \details There are no z - derivatives in a 2D case

template < class DataType, int DIM >
void FELagrangeQuad< DataType, DIM >::N_yz(const Coord &pt,
                                      std::vector< DataType > &weight) const {
  assert (DIM == 2);
  weight.assign(weight.size(), 0.);
}

/// \details There are no z - derivatives in a 2D case

template < class DataType, int DIM >
void FELagrangeQuad< DataType, DIM >::N_zz(const Coord &pt,
                                      std::vector< DataType > &weight) const {
  assert (DIM == 2);
  weight.assign(weight.size(), 0.);
}

} // namespace doffem
} // namespace hiflow
#endif
