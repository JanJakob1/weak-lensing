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

#ifndef __FEM_FELAGRANGE_HEX_H_
#define __FEM_FELAGRANGE_HEX_H_

#include "felagrange.h"
#include <cassert>

namespace hiflow {
namespace doffem {

///
/// \class FELagrangeHex felagrange_hex.h
/// \brief Lagrangian Finite Element on a Hexahedron
/// \author Michael Schick<br>Martin Baumann
///

template < class DataType, int DIM >
class FELagrangeHex : public virtual FELagrange< DataType, DIM > {
public:
  typedef Vec<DIM, DataType> Coord;

  /// Default Constructor
  FELagrangeHex();

  /// Default Destructor
  virtual ~FELagrangeHex();

  std::string get_name() const { return "LagrangeHexahedron"; }

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
  /// Index ordering form space (x=i,y=j,z=k) to vector index
  int ijk2ind(int i, int j, int k) const;
};

template < class DataType, int DIM > FELagrangeHex< DataType, DIM >::FELagrangeHex() {
  this->my_type_ = FEType< DataType, DIM >::LAGRANGE_HEX;

  // initialize reference cell

  assert(this->ref_cell_ == NULL);
  this->ref_cell_ = &(mesh::CellType::get_instance(mesh::CellType::HEXAHEDRON));
}

template < class DataType, int DIM > FELagrangeHex< DataType, DIM >::~FELagrangeHex() {}

template < class DataType, int DIM > void FELagrangeHex< DataType, DIM >::init_coord() {
  assert (DIM == 3);
  
  // set topological degree

  this->tdim_ = 3;

  // set coordinates of all DoF points in lexicographical ordering

  if (this->my_deg_ == 0) {
    this->coord_.clear();

    // Coordinates of middle-point of hex
    Coord coord;

    coord[0] = 0.5;
    coord[1] = 0.5;
    coord[2] = 0.5;

    this->coord_.push_back(coord);
  } else {
    assert(this->my_deg_ > 0);

    const DataType offset = (1.0 / this->my_deg_);
    const int fd_1 = this->my_deg_ + 1;

    const int nb_dof_on_cell = fd_1 * fd_1 * fd_1;
    const int nb_dof_on_line = fd_1;

    this->coord_.clear();
    this->coord_.resize(nb_dof_on_cell);

    /*     Triquadratic Hexaedron (Q2-Element)
           ===================================

                24----------25-----------26
                /|          /|          /|
               / |         / |         / |
             15----------16-----------17 |
             /|  |       /|  |       /|  |
            / | 21------/---22------/-|--23
           6-----------7-----------8  | /|
           |  |/ |     |  |/ |     |  |/ |
           | 12--|-----|-13--|-----|--14 |
           | /|  |     | /|  |     | /|  |
           |/ | 18-----|----19-----|/-|--20
           3-----------4-----------5  | /
           |  |/       |  |/       |  |/
           |  9--------|-10--------|--11
           | /         | /         | /
           |/          |/          |/
           0-----------1-----------2          */

    // Filling coord vector for full hex by lexicographical strategy

    for (int k = 0; k < nb_dof_on_line; ++k) {
      const DataType k_offset = k * offset;
      for (int j = 0; j < nb_dof_on_line; ++j) {
        const DataType j_offset = j * offset;
        for (int i = 0; i < nb_dof_on_line; ++i) {
          Coord coord;

          coord[0] = i * offset;
          coord[1] = j_offset;
          coord[2] = k_offset;

          this->coord_[ijk2ind(i, j, k)] = coord;
        }
      }
    }
  }
}

template < class DataType, int DIM >
int FELagrangeHex< DataType, DIM >::ijk2ind(int i, int j, int k) const {
  assert (DIM == 3);
  const int nb_dof_on_line = this->my_deg_ + 1;
  return (i + (j + k * nb_dof_on_line) * nb_dof_on_line);
}

/// \details Every degree of a used lagrangian polynomial has to satisfy the
///          condition degree < my_deg_. All possible combinations are
///          multiplied and considered in the sum.

template < class DataType, int DIM >
void FELagrangeHex< DataType, DIM >::N(const Coord &pt,
                                  std::vector< DataType > &weight) const {
  assert (DIM == 3);
  assert(weight.size() == this->get_nb_dof_on_cell());

  if (this->my_deg_ > 0) {
    for (int k = 0; k <= this->my_deg_; ++k) {
      const DataType lp_k = this->lp_.poly(this->my_deg_, k, pt[2]);
      for (int j = 0; j <= this->my_deg_; ++j) {
        const DataType lp_j = this->lp_.poly(this->my_deg_, j, pt[1]);
        for (int i = 0; i <= this->my_deg_; ++i) {
          weight[ijk2ind(i, j, k)] =
              this->lp_.poly(this->my_deg_, i, pt[0]) * lp_j * lp_k;
        }
      }
    }
  } else
    weight[ijk2ind(0, 0, 0)] = 1.0;
}

/// \details Every degree of a used lagrangian polynomial has to satisfy the
///          condition degree < my_deg_. All possible combinations are
///          multiplied and considered in the sum, w.r.t. the derivatives for
///          the x - variable.

template < class DataType, int DIM >
void FELagrangeHex< DataType, DIM >::N_x(const Coord &pt,
                                    std::vector< DataType > &weight) const {
  assert (DIM == 3);
  assert(weight.size() == this->get_nb_dof_on_cell());

  if (this->my_deg_ > 0) {
    for (int k = 0; k <= this->my_deg_; ++k) {
      const DataType lp_k = this->lp_.poly(this->my_deg_, k, pt[2]);
      for (int j = 0; j <= this->my_deg_; ++j) {
        const DataType lp_j = this->lp_.poly(this->my_deg_, j, pt[1]);
        for (int i = 0; i <= this->my_deg_; ++i) {
          weight[ijk2ind(i, j, k)] =
              this->lp_.poly_x(this->my_deg_, i, pt[0]) * lp_j * lp_k;
        }
      }
    }
  } else {
    weight[ijk2ind(0, 0, 0)] = 0.0;
  }
}

/// \details Every degree of a used lagrangian polynomial has to satisfy the
///          condition degree < my_deg_. All possible combinations are
///          multiplied and considered in the sum, w.r.t. the derivatives for
///          the y - variable.

template < class DataType, int DIM >
void FELagrangeHex< DataType, DIM >::N_y(const Coord &pt,
                                    std::vector< DataType > &weight) const {
  assert (DIM == 3);
  assert(weight.size() == this->get_nb_dof_on_cell());

  if (this->my_deg_ > 0) {
    for (int k = 0; k <= this->my_deg_; ++k) {
      const DataType lp_k = this->lp_.poly(this->my_deg_, k, pt[2]);
      for (int j = 0; j <= this->my_deg_; ++j) {
        const DataType lp_j = this->lp_.poly_x(this->my_deg_, j, pt[1]);
        for (int i = 0; i <= this->my_deg_; ++i) {
          weight[ijk2ind(i, j, k)] =
              this->lp_.poly(this->my_deg_, i, pt[0]) * lp_j * lp_k;
        }
      }
    }
  } else {
    weight[ijk2ind(0, 0, 0)] = 0.0;
  }
}

/// \details Every degree of a used lagrangian polynomial has to satisfy the
///          condition degree < my_deg_. All possible combinations are
///          multiplied and considered in the sum, w.r.t. the derivatives for
///          the z - variable.

template < class DataType, int DIM >
void FELagrangeHex< DataType, DIM >::N_z(const Coord &pt,
                                    std::vector< DataType > &weight) const {
  assert (DIM == 3);
  assert(weight.size() == this->get_nb_dof_on_cell());

  if (this->my_deg_ > 0) {
    for (int k = 0; k <= this->my_deg_; ++k) {
      const DataType lp_k = this->lp_.poly_x(this->my_deg_, k, pt[2]);
      for (int j = 0; j <= this->my_deg_; ++j) {
        const DataType lp_j = this->lp_.poly(this->my_deg_, j, pt[1]);
        for (int i = 0; i <= this->my_deg_; ++i) {
          weight[ijk2ind(i, j, k)] =
              this->lp_.poly(this->my_deg_, i, pt[0]) * lp_j * lp_k;
        }
      }
    }
  } else {
    weight[ijk2ind(0, 0, 0)] = 0.0;
  }
}

/// \details Every degree of a used lagrangian polynomial has to satisfy the
///          condition degree < my_deg_. All possible combinations are
///          multiplied and considered in the sum, w.r.t. the second derivatives
///          for the xx - variable.

template < class DataType, int DIM >
void FELagrangeHex< DataType, DIM >::N_xx(const Coord &pt,
                                     std::vector< DataType > &weight) const {
  assert (DIM == 3);
  assert(weight.size() == this->get_nb_dof_on_cell());

  if (this->my_deg_ > 0) {
    for (int k = 0; k <= this->my_deg_; ++k) {
      const DataType lp_k = this->lp_.poly(this->my_deg_, k, pt[2]);
      for (int j = 0; j <= this->my_deg_; ++j) {
        const DataType lp_j = this->lp_.poly(this->my_deg_, j, pt[1]);
        for (int i = 0; i <= this->my_deg_; ++i) {
          weight[ijk2ind(i, j, k)] =
              this->lp_.poly_xx(this->my_deg_, i, pt[0]) * lp_j * lp_k;
        }
      }
    }
  } else {
    weight[ijk2ind(0, 0, 0)] = 0.0;
  }
}

/// \details Every degree of a used lagrangian polynomial has to satisfy the
///          condition degree < my_deg_. All possible combinations are
///          multiplied and considered in the sum, w.r.t. the second derivatives
///          for the xy - variable.

template < class DataType, int DIM >
void FELagrangeHex< DataType, DIM >::N_xy(const Coord &pt,
                                     std::vector< DataType > &weight) const {
  assert (DIM == 3);
  assert(weight.size() == this->get_nb_dof_on_cell());

  if (this->my_deg_ > 0) {
    for (int k = 0; k <= this->my_deg_; ++k) {
      const DataType lp_k = this->lp_.poly(this->my_deg_, k, pt[2]);
      for (int j = 0; j <= this->my_deg_; ++j) {
        const DataType lp_j = this->lp_.poly_x(this->my_deg_, j, pt[1]);
        for (int i = 0; i <= this->my_deg_; ++i) {
          weight[ijk2ind(i, j, k)] =
              this->lp_.poly_x(this->my_deg_, i, pt[0]) * lp_j * lp_k;
        }
      }
    }
  } else {
    weight[ijk2ind(0, 0, 0)] = 0.0;
  }
}

/// \details Every degree of a used lagrangian polynomial has to satisfy the
///          condition degree < my_deg_. All possible combinations are
///          multiplied and considered in the sum, w.r.t. the second derivatives
///          for the xz - variable.

template < class DataType, int DIM >
void FELagrangeHex< DataType, DIM >::N_xz(const Coord &pt,
                                     std::vector< DataType > &weight) const {
  assert (DIM == 3);
  assert(weight.size() == this->get_nb_dof_on_cell());

  if (this->my_deg_ > 0) {
    for (int k = 0; k <= this->my_deg_; ++k) {
      const DataType lp_k = this->lp_.poly_x(this->my_deg_, k, pt[2]);
      for (int j = 0; j <= this->my_deg_; ++j) {
        const DataType lp_j = this->lp_.poly(this->my_deg_, j, pt[1]);
        for (int i = 0; i <= this->my_deg_; ++i) {
          weight[ijk2ind(i, j, k)] =
              this->lp_.poly_x(this->my_deg_, i, pt[0]) * lp_j * lp_k;
        }
      }
    }
  } else {
    weight[ijk2ind(0, 0, 0)] = 0.0;
  }
}

/// \details Every degree of a used lagrangian polynomial has to satisfy the
///          condition degree < my_deg_. All possible combinations are
///          multiplied and considered in the sum, w.r.t. the second derivatives
///          for the yy - variable.

template < class DataType, int DIM >
void FELagrangeHex< DataType, DIM >::N_yy(const Coord &pt,
                                     std::vector< DataType > &weight) const {
  assert (DIM == 3);
  assert(weight.size() == this->get_nb_dof_on_cell());

  if (this->my_deg_ > 0) {
    for (int k = 0; k <= this->my_deg_; ++k) {
      const DataType lp_k = this->lp_.poly(this->my_deg_, k, pt[2]);
      for (int j = 0; j <= this->my_deg_; ++j) {
        const DataType lp_j = this->lp_.poly_xx(this->my_deg_, j, pt[1]);
        for (int i = 0; i <= this->my_deg_; ++i) {
          weight[ijk2ind(i, j, k)] =
              this->lp_.poly(this->my_deg_, i, pt[0]) * lp_j * lp_k;
        }
      }
    }
  } else {
    weight[ijk2ind(0, 0, 0)] = 0.0;
  }
}

/// \details Every degree of a used lagrangian polynomial has to satisfy the
///          condition degree < my_deg_. All possible combinations are
///          multiplied and considered in the sum, w.r.t. the second derivatives
///          for the yz - variable.

template < class DataType, int DIM >
void FELagrangeHex< DataType, DIM >::N_yz(const Coord &pt,
                                     std::vector< DataType > &weight) const {
  assert (DIM == 3);
  assert(weight.size() == this->get_nb_dof_on_cell());

  if (this->my_deg_ > 0) {
    for (int k = 0; k <= this->my_deg_; ++k) {
      const DataType lp_k = this->lp_.poly_x(this->my_deg_, k, pt[2]);
      for (int j = 0; j <= this->my_deg_; ++j) {
        const DataType lp_j = this->lp_.poly_x(this->my_deg_, j, pt[1]);
        for (int i = 0; i <= this->my_deg_; ++i) {
          weight[ijk2ind(i, j, k)] =
              this->lp_.poly(this->my_deg_, i, pt[0]) * lp_j * lp_k;
        }
      }
    }
  } else {
    weight[ijk2ind(0, 0, 0)] = 0.0;
  }
}

/// \details Every degree of a used lagrangian polynomial has to satisfy the
///          condition degree < my_deg_. All possible combinations are
///          multiplied and considered in the sum, w.r.t. the second derivatives
///          for the zz - variable.

template < class DataType, int DIM >
void FELagrangeHex< DataType, DIM >::N_zz(const Coord &pt,
                                     std::vector< DataType > &weight) const {
  assert (DIM == 3);
  assert(weight.size() == this->get_nb_dof_on_cell());

  if (this->my_deg_ > 0) {
    for (int k = 0; k <= this->my_deg_; ++k) {
      const DataType lp_k = this->lp_.poly_xx(this->my_deg_, k, pt[2]);
      for (int j = 0; j <= this->my_deg_; ++j) {
        const DataType lp_j = this->lp_.poly(this->my_deg_, j, pt[1]);
        for (int i = 0; i <= this->my_deg_; ++i) {
          weight[ijk2ind(i, j, k)] =
              this->lp_.poly(this->my_deg_, i, pt[0]) * lp_j * lp_k;
        }
      }
    }
  } else {
    weight[ijk2ind(0, 0, 0)] = 0.0;
  }
}

} // namespace doffem
} // namespace hiflow
#endif
