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

#ifndef _FEM_FELAGRANGE_H_
#define _FEM_FELAGRANGE_H_

#include "fem/fetype.h"
#include "polynomials/lagrangepolynomial.h"

namespace hiflow {
namespace doffem {

///
/// \class FELagrange felagrange.h
/// \brief Lagrangian Finite Element
/// \author Michael Schick<br>Martin Baumann
///

template < class DataType, int DIM > 
class FELagrange : public FEType< DataType, DIM > {
public:
  typedef Vec<DIM, DataType> Coord;

  /// Default Constructor
  FELagrange();

  /// For given point, get values of all shapefunctions on reference cell
  virtual void N(const Coord &pt, std::vector< DataType > &weight) const = 0;
  void N(const Coord &pt, int comp, std::vector< DataType > &weight) const {
    assert (comp == 0);
    this->N(pt, weight);
  }

  /// For given coordinates, get x - derivative of all shapefunctions on
  /// reference cell
  virtual void N_x(const Coord &pt, std::vector< DataType > &weight) const = 0;
  void N_x(const Coord &pt, int comp, std::vector< DataType > &weight) const {
    assert (comp == 0);
    this->N_x(pt, weight);
  }
  
  /// For given coordinates, get y - derivative of all shapefunctions on
  /// reference cell
  virtual void N_y(const Coord &pt, std::vector< DataType > &weight) const = 0;
  void N_y(const Coord &pt, int comp, std::vector< DataType > &weight) const {
    assert (comp == 0);
    this->N_y(pt, weight);
  }

  /// For given coordinates, get z - derivative of all shapefunctions on
  /// reference cell
  virtual void N_z(const Coord &pt, std::vector< DataType > &weight) const = 0;
  void N_z(const Coord &pt, int comp, std::vector< DataType > &weight) const {
    assert (comp == 0);
    this->N_z(pt, weight);
  }

  /// For given coordinates, get xx - derivative of all shapefunctions on
  /// reference cell
  virtual void N_xx(const Coord &pt, std::vector< DataType > &weight) const = 0;
  void N_xx(const Coord &pt, int comp, std::vector< DataType > &weight) const {
    assert (comp == 0);
    this->N_xx(pt, weight);
  }

  /// For given coordinates, get xy - derivative of all shapefunctions on
  /// reference cell
  virtual void N_xy(const Coord &pt, std::vector< DataType > &weight) const = 0;
  void N_xy(const Coord &pt, int comp, std::vector< DataType > &weight) const {
    assert (comp == 0);
    this->N_xy(pt, weight);
  }

  /// For given coordinates, get xz - derivative of all shapefunctions on
  /// reference cell
  virtual void N_xz(const Coord &pt, std::vector< DataType > &weight) const = 0;
  void N_xz(const Coord &pt, int comp, std::vector< DataType > &weight) const {
    assert (comp == 0);
    this->N_xz(pt, weight);
  }

  /// For given coordinates, get yy - derivative of all shapefunctions on
  /// reference cell
  virtual void N_yy(const Coord &pt, std::vector< DataType > &weight) const = 0;
  void N_yy(const Coord &pt, int comp, std::vector< DataType > &weight) const {
    assert (comp == 0);
    this->N_yy(pt, weight);
  }

  /// For given coordinates, get yz - derivative of all shapefunctions on
  /// reference cell
  virtual void N_yz(const Coord &pt, std::vector< DataType > &weight) const = 0;
  void N_yz(const Coord &pt, int comp, std::vector< DataType > &weight) const {
    assert (comp == 0);
    this->N_yz(pt, weight);
  }

  /// For given coordinates, get zz - derivative of all shapefunctions on
  /// reference cell
  virtual void N_zz(const Coord &pt, std::vector< DataType > &weight) const = 0;
  void N_zz(const Coord &pt, int comp, std::vector< DataType > &weight) const {
    assert (comp == 0);
    this->N_zz(pt, weight);
  }

protected:
  /// Lagrange polynomials which are used for evaluating shapefunctions
  LagrangePolynomial< DataType > lp_;
};

template < class DataType, int DIM > 
FELagrange< DataType, DIM >::FELagrange() {
  this->nb_comp_ = 1;
}

} // namespace doffem
} // namespace hiflow
#endif
