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

#ifndef __FEM_CELL_TRANSFORMATION_H_
#define __FEM_CELL_TRANSFORMATION_H_

#include <cassert>
#include <vector>
#include <limits>
#include <cmath>

#include "common/log.h"
#include "common/macros.h"
#include "common/vector_algebra.h"
#include "dof/dof_fem_types.h"
#include "fem/reference_cell.h"

namespace hiflow {
namespace doffem {

enum CellTrafoType
{
  CELL_TRAFO_NONE = -1,
  CELL_TRAFO_LINEARLINE = 0,
  CELL_TRAFO_LINEARTRI = 1,
  CELL_TRAFO_LINEARTET = 2,
  CELL_TRAFO_LINEARPYR = 3,
  CELL_TRAFO_BILINEARQUAD = 4,
  CELL_TRAFO_TRILINEARHEX = 5,
  CELL_TRAFO_ALIGNEDQUAD = 6,
  CELL_TRAFO_ALIGNEDHEX = 7
};

///
/// \class CellTransformation cell_transformation.h
/// \brief Ancestor class of all transformation mappings from reference to
/// physical cells 
/// \author Michael Schick<br>Martin Baumann<br>Simon Gawlok<br>Philipp Gerstner
///

template < class DataType, int DIM > 
class CellTransformation 
{
public:
  typedef Vec<DIM, DataType> Coord;

  /// Use this constructor which needs the geometrical dimension as input
  explicit CellTransformation(ConstRefCellPtr<DataType, DIM> ref_cell);

  virtual ~CellTransformation() 
  {
  }

  std::string name() const 
  {
    return this->name_;
  }
  
#if 0
  virtual void copy_from (CellTransformation<DataType, DIM> const * rhs)
  {
    // Copy only works for celltrafos os same type
    assert (this->trafo_type_ == rhs->trafo_type_);
    assert (rhs != nullptr);
    this->ref_cell_ = rhs->ref_cell_;
    this->fixed_ref_cell_type_ = rhs->fixed_ref_cell_type_; 
    this->coord_vtx_ = rhs->coord_vtx_;
    this->order_ = rhs->order_;
  } 
#endif

  virtual bool differs_by_translation_from (ConstCellTrafoPtr<DataType, DIM> rhs) const = 0;
  
  /// Reinitialization of the transformation via coordinates of physical cell
  virtual void reinit(const std::vector<DataType> &coord_vtx);

  ConstRefCellPtr<DataType, DIM> get_ref_cell() const
  {
    return this->ref_cell_;
  }
  
  std::vector< Coord > get_reference_coordinates() const 
  { 
    assert(this->ref_cell_);
    return this->ref_cell_->get_coords();
  }

  std::vector<Coord> get_coordinates() const 
  { 
    return this->coord_vtx_;
  }

  void print_vertex_coords() const 
  {
    for (size_t l=0; l<this->coord_vtx_.size(); ++l)
    {
      std::cout << "[" << l << "]: " << this->coord_vtx_[l] << std::endl;
    }
    std::cout << std::endl;
  }
  
  /// \brief Check whether a given point is contained in the closure of the
  /// cell.
  ///
  /// \param[in] coord_ref    reference coordinates of the point
  /// \returns  True if reference coordinates are contained in the cell.
  bool contains_reference_point(const Coord &coord_ref) const
  {
    assert(this->ref_cell_);
    return this->ref_cell_->contains_point (coord_ref);
  }

  /// \brief Check whether a given point is contained in the closure of the
  /// cell.
  /// \param[in]  coord_phys   Physical coordinates of the point.
  /// \param[out] coord_ref    Optional output parameter, where reference
  /// coordinates in the cell are stored, if the function returns true. \returns
  /// True if reference coordinates are contained in the cell.
  virtual bool contains_physical_point(const Coord &coord_phys, Coord &coord_ref) const;

  void transform ( const Coord& coord_ref, Coord& coord_mapped ) const;
  
  Coord transform ( const Coord& coord_ref ) const
  {
    Coord coord_mapped;
    this->transform(coord_ref, coord_mapped);
    return coord_mapped;
  }

  void J (const Coord &coord_ref, Mat<DIM, DIM, DataType>& J) const;

  void H (const Coord &coord_ref, size_t d, Mat<DIM, DIM, DataType>& H) const;

  /// \brief compute determinant of Jacobian of cell transformation at given reference coordinates
  DataType detJ (const Coord &coord_ref) const;

  void J_and_detJ (const Coord &coord_ref, Mat<DIM, DIM, DataType>& J, DataType& detJ) const;
  
  /// \brief compute gradient of determinant of Jacobian of cell transformation at given reference coordinates
  void grad_detJ (const Coord &coord_ref, Vec<DIM, DataType> & grad) const;
  void grad_inv_detJ (const Coord &coord_ref, Vec<DIM, DataType> & grad) const;
  
  DataType detJ_x (const Coord &coord_ref) const;
  DataType detJ_y (const Coord &coord_ref) const;
  DataType detJ_z (const Coord &coord_ref) const;
  
  /// \brief compute hessian ofdeterminant of Jacobian of cell transformation at given reference coordinates
  void hessian_detJ (const Coord &coord_ref, Mat<DIM, DIM, DataType> & mat) const;
  DataType detJ_xx (const Coord &coord_ref) const;
  DataType detJ_xy (const Coord &coord_ref) const;
  DataType detJ_xz (const Coord &coord_ref) const;
  DataType detJ_yy (const Coord &coord_ref) const;
  DataType detJ_yz (const Coord &coord_ref) const;          
  DataType detJ_zz (const Coord &coord_ref) const;
  
  /// \brief Given physical cell coordinates in 1D,
  ///        this routine computes the corresponding reference cell coordinates
  /// @return true, if inverse computation was successful
  virtual bool inverse(Coord co_phy, Coord &co_ref) const = 0;

  /// Given reference coordinates, this routine computes the physical x
  /// coordinates
  virtual DataType x(const Coord &coord_ref) const { return 0.; }
  
  /// \brief Given reference coordinates, this routine computes the derivatives
  /// of the mapping (ref_coordinates to physical x value)
  virtual DataType x_x(const Coord &coord_ref) const { return 0.; }
  virtual DataType x_y(const Coord &coord_ref) const { return 0.; }
  virtual DataType x_z(const Coord &coord_ref) const { return 0.; }
  
  /// \brief Given reference coordinates, these routine compute the second
  /// derivatives of the mapping (ref_coordinates to physical x value). 
  virtual DataType x_xx(const Coord &coord_ref) const { return 0.; }
  virtual DataType x_xy(const Coord &coord_ref) const { return 0.; }
  virtual DataType x_xz(const Coord &coord_ref) const { return 0.; }
  virtual DataType x_yy(const Coord &coord_ref) const { return 0.; }
  virtual DataType x_yz(const Coord &coord_ref) const { return 0.; }
  virtual DataType x_zz(const Coord &coord_ref) const { return 0.; }

  /// \brief Given reference coordinates, these routine compute the third
  /// derivatives of the mapping (ref_coordinates to physical x value). 
  virtual DataType x_xxx(const Coord &coord_ref) const { return 0.; }
  virtual DataType x_xxy(const Coord &coord_ref) const { return 0.; }
  virtual DataType x_xxz(const Coord &coord_ref) const { return 0.; }
  virtual DataType x_xyy(const Coord &coord_ref) const { return 0.; }
  virtual DataType x_xyz(const Coord &coord_ref) const { return 0.; }
  virtual DataType x_xzz(const Coord &coord_ref) const { return 0.; }
  virtual DataType x_yyy(const Coord &coord_ref) const { return 0.; }
  virtual DataType x_yyz(const Coord &coord_ref) const { return 0.; }
  virtual DataType x_yzz(const Coord &coord_ref) const { return 0.; }
  virtual DataType x_zzz(const Coord &coord_ref) const { return 0.; }
  
  /// Given reference coordinates, this computes the physical y coordinates
  virtual DataType y(const Coord &coord_ref) const { return 0.; }
  
  /// \brief Given reference coordinates, this routine computes the derivatives
  /// of the mapping (ref_coordinates to physical y value)
  virtual DataType y_x(const Coord &coord_ref) const { return 0.; }
  virtual DataType y_y(const Coord &coord_ref) const { return 0.; }
  virtual DataType y_z(const Coord &coord_ref) const { return 0.; }

  /// \brief Given reference coordinates, this routine computes the second
  /// derivatives of the mapping (ref_coordinates to physical y value)
  virtual DataType y_xx(const Coord &coord_ref) const { return 0.; }
  virtual DataType y_xy(const Coord &coord_ref) const { return 0.; }
  virtual DataType y_xz(const Coord &coord_ref) const { return 0.; }
  virtual DataType y_yy(const Coord &coord_ref) const { return 0.; }
  virtual DataType y_yz(const Coord &coord_ref) const { return 0.; }
  virtual DataType y_zz(const Coord &coord_ref) const { return 0.; }

  /// \brief Given reference coordinates, these routine compute the third
  /// derivatives of the mapping (ref_coordinates to physical y value). 
  virtual DataType y_xxx(const Coord &coord_ref) const { return 0.; }
  virtual DataType y_xxy(const Coord &coord_ref) const { return 0.; }
  virtual DataType y_xxz(const Coord &coord_ref) const { return 0.; }
  virtual DataType y_xyy(const Coord &coord_ref) const { return 0.; }
  virtual DataType y_xyz(const Coord &coord_ref) const { return 0.; }
  virtual DataType y_xzz(const Coord &coord_ref) const { return 0.; }
  virtual DataType y_yyy(const Coord &coord_ref) const { return 0.; }
  virtual DataType y_yyz(const Coord &coord_ref) const { return 0.; }
  virtual DataType y_yzz(const Coord &coord_ref) const { return 0.; }
  virtual DataType y_zzz(const Coord &coord_ref) const { return 0.; }

  /// Given reference coordinates, this computes the physical z coordinates
  virtual DataType z(const Coord &coord_ref) const { return 0.; }
  /// \brief Given reference coordinates, this routine computes the derivatives
  /// in of the mapping (ref_coordinates to physical z value).
  /// Return value is a dummy value for 2D problems, but neccessary for UnitIntegrator

  /// \brief Given reference coordinates, this routine computes the derivatives
  /// of the mapping (ref_coordinates to physical z value).
  /// Return value is a dummy value for 2D problems, but neccessary for
  /// UnitIntegrator
  virtual DataType z_x(const Coord &coord_ref) const { return 0.; }
  virtual DataType z_y(const Coord &coord_ref) const { return 0.; }
  virtual DataType z_z(const Coord &coord_ref) const { return 0.; }

  /// \brief Given reference coordinates, this routine computes the second
  /// derivatives  of the mapping (ref_coordinates to physical z
  /// value). Return value is a dummy value for 2D problems, but
  /// neccessary for UnitIntegrator
  virtual DataType z_xx(const Coord &coord_ref) const { return 0.; }
  virtual DataType z_xy(const Coord &coord_ref) const { return 0.; }
  virtual DataType z_xz(const Coord &coord_ref) const { return 0.; }
  virtual DataType z_yy(const Coord &coord_ref) const { return 0.; }
  virtual DataType z_yz(const Coord &coord_ref) const { return 0.; }
  virtual DataType z_zz(const Coord &coord_ref) const { return 0.; }

  /// \brief Given reference coordinates, these routine compute the third
  /// derivatives of the mapping (ref_coordinates to physical x value). 
  virtual DataType z_xxx(const Coord &coord_ref) const { return 0.; }
  virtual DataType z_xxy(const Coord &coord_ref) const { return 0.; }
  virtual DataType z_xxz(const Coord &coord_ref) const { return 0.; }
  virtual DataType z_xyy(const Coord &coord_ref) const { return 0.; }
  virtual DataType z_xyz(const Coord &coord_ref) const { return 0.; }
  virtual DataType z_xzz(const Coord &coord_ref) const { return 0.; }
  virtual DataType z_yyy(const Coord &coord_ref) const { return 0.; }
  virtual DataType z_yyz(const Coord &coord_ref) const { return 0.; }
  virtual DataType z_yzz(const Coord &coord_ref) const { return 0.; }
  virtual DataType z_zzz(const Coord &coord_ref) const { return 0.; }

protected:
  /// \details The vector index is calculated by an offset of the
  ///          magnitude i * geometrical dimension
  /// \param[in] i index of vertex id
  /// \param[in] j index of coordinate id (0 for x, 1 for y and 2 for z)
  /// \return index for vector of coordinates coord_vtx_

  inline int ij2ind(int i, int j) const 
  { 
    return i * DIM + j; 
  }

  bool inverse_newton( Coord &co_phy, Coord &co_ref ) const;
  
  bool inverse_newton( Coord &co_phy, Coord &co_ref, Coord &co_ref_0 ) const;
  
  bool inverse_newton_2Dto3D( Vec<3, DataType> &co_phy, Vec<2, DataType> &co_ref, Vec<2, DataType> &co_ref_0) const;

  ConstRefCellPtr<DataType, DIM> ref_cell_;
  RefCellType fixed_ref_cell_type_;

  CellTrafoType trafo_type_;
  
  /// Vector, which holds the coordinates of every vertex of the physical cell
  std::vector<Coord> coord_vtx_;
  
  /// highest polynomial order 
  int order_;
  
  std::string name_;
};

template < class DataType, int DIM >
CellTransformation< DataType, DIM >::CellTransformation(ConstRefCellPtr<DataType, DIM> ref_cell) 
: ref_cell_(ref_cell), 
order_(10),
trafo_type_(CELL_TRAFO_NONE) 
{
}

/// \details Given vector of coordinates on physical cell, the are
///          stored the protected member variable coord_vtx_

template < class DataType, int DIM >
void CellTransformation< DataType, DIM >::reinit(const std::vector<DataType> &coord_vtx) 
{  
  assert (this->fixed_ref_cell_type_ == this->ref_cell_->type());
  
  int num_points = coord_vtx.size() / DIM;
  coord_vtx_.clear();
  coord_vtx_.resize(num_points);
  for (int i=0; i<num_points; ++i) 
  {
    for (int d=0; d<DIM; ++d) 
    {
      coord_vtx_[i][d] = coord_vtx[i*DIM + d];
    } 
  }
}

/// \details The Inverse of the mapping (reference cell to physical cell)
///          is computed in 3D via a Newton-scheme with Armijo updates. The
///          Problem, which is to be solved writes<br>
///          G([x_ref, y_ref, z_ref]) := [x_phy - x(coord_ref),
///                                       y_phy - y(coord_ref),
///                                       z_phy - z(coord_ref)]<br>
///          solve G([x_ref, y_ref, z_ref]) = 0 via Newton-scheme
/// \param[in] x_phy x coordinate on physical cell
/// \param[in] y_phy y coordinate on physical cell
/// \param[in] z_phy z coordinate on physical cell
/// \param[out] x_ref x coordinate on reference cell
/// \param[out] y_ref y coordinate on reference cell
/// \param[out] z_ref z coordinate on reference cell
template < class DataType, int DIM >
bool CellTransformation< DataType, DIM >::inverse_newton( Coord &co_phy, Coord &co_ref ) const
{
  Coord co_ref_0;
  switch (DIM)
  {
    case 1:
      co_ref_0[0] = 0.5;
      break;
    case 2:
      co_ref_0[0] = 0.55;
      co_ref_0[1] = 0.55;
      break;
    case 3:
      co_ref_0[0] = 0.1154;
      co_ref_0[1] = 0.1832;
      co_ref_0[2] = 0.1385;
      break;
  }
  return this->inverse_newton (co_phy, co_ref, co_ref_0);
}

template < class DataType, int DIM >
bool CellTransformation< DataType, DIM >::inverse_newton( Coord &co_phy, Coord &co_ref, Coord &co_ref_0) const {
  // Initialisation

  Coord pt_phy = co_phy;

  Coord ref_k1, ref_k;
  ref_k = co_ref_0;

  const DataType tol_eps = 1.e3 * std::numeric_limits< DataType >::epsilon();

  // Some general parameters

  const int iter_max = 1000;
  int iter = 0;

  // Residual

  Coord coord_ref_k = ref_k;
  Coord pt_k;
  switch (DIM)
  {
    case 1:
      pt_k[0] = this->x(coord_ref_k);
      break;
    case 2:
      pt_k[0] = this->x(coord_ref_k);
      pt_k[1] = this->y(coord_ref_k);
      break;
    case 3:
      pt_k[0] = this->x(coord_ref_k);
      pt_k[1] = this->y(coord_ref_k);
      pt_k[2] = this->z(coord_ref_k);
      break;
  }
  
  DataType residual = norm(pt_phy - pt_k);
  DataType abserr = norm(ref_k1 - ref_k);

  // Newton

  // Jacobian Matrix (grad G)
  Mat< DIM, DIM, DataType > G;
  // Inverse of the jacobian
  Mat< DIM, DIM, DataType > B;

  while (residual > 10. * std::numeric_limits< DataType >::epsilon()) 
  {
    switch (DIM)
    {
      case 1:
        G(0, 0) = this->x_x(coord_ref_k);
        break;
      case 2:
        G(0, 0) = this->x_x(coord_ref_k);
        G(0, 1) = this->x_y(coord_ref_k);
    
        G(1, 0) = this->y_x(coord_ref_k);
        G(1, 1) = this->y_y(coord_ref_k);
        break;
      case 3:
        G(0, 0) = this->x_x(coord_ref_k);
        G(0, 1) = this->x_y(coord_ref_k);
        G(0, 2) = this->x_z(coord_ref_k);

        G(1, 0) = this->y_x(coord_ref_k);
        G(1, 1) = this->y_y(coord_ref_k);
        G(1, 2) = this->y_z(coord_ref_k);

        G(2, 0) = this->z_x(coord_ref_k);
        G(2, 1) = this->z_y(coord_ref_k);
        G(2, 2) = this->z_z(coord_ref_k);
    }
  
#ifndef NDEBUG
    const DataType detG = det(G);

    assert(detG != 0.);
#endif

    inv(G, B);

    // Armijo parameter

    const int iter_max_armijo = 500;
    int iter_armijo = 0;

    DataType residual_armijo = 2. * residual;

    DataType omega = 1.;

    const Coord update_vec = B * (pt_phy - pt_k);

    // Start Armijo

    while ((iter_armijo <= iter_max_armijo) && (residual_armijo > residual)) {

      ref_k1 = ref_k + omega * update_vec;

      for (int d=0; d<DIM; ++d) 
      {
        if ((ref_k1[d] >= -tol_eps) && (ref_k1[d] <= tol_eps)) {
          ref_k1[d] = 0.;
        }
        if ((ref_k1[d] - 1. >= -tol_eps) && (ref_k1[d] - 1. <= tol_eps)) {
          ref_k1[d] = 1.;
        }
      }

      Coord ref_check = ref_k1;

      while (!(this->contains_reference_point(ref_check))) {
        omega /= 2.0;
        ref_k1 = ref_k + omega * update_vec;

        for (int d=0; d<DIM; ++d) 
        {
          if ((ref_k1[d] >= -tol_eps) && (ref_k1[d] <= tol_eps)) {
            ref_k1[d] = 0.;
          }
          if ((ref_k1[d] - 1. >= -tol_eps) && (ref_k1[d] - 1. <= tol_eps)) {
            ref_k1[d] = 1.;
          }
        }
        ref_check = ref_k1;
      }

      Coord coord_ref_k1 = ref_k1;

      Coord F_k1;
      switch (DIM)
      {
        case 1:
          F_k1[0] = this->x(coord_ref_k1);
          break;
        case 2:
          F_k1[0] = this->x(coord_ref_k1);
          F_k1[1] = this->y(coord_ref_k1);
          break;
        case 3:
          F_k1[0] = this->x(coord_ref_k1);
          F_k1[1] = this->y(coord_ref_k1);
          F_k1[2] = this->z(coord_ref_k1);
          break;
      }

      residual_armijo = norm(pt_phy - F_k1);

      ++iter_armijo;
      omega /= 2.;
    }

    abserr = norm(ref_k1 - ref_k);
    ref_k = ref_k1;

    coord_ref_k = ref_k;

    switch (DIM)
    {
      case 1:
        pt_k[0] = this->x(coord_ref_k);
        break;
      case 2:
        pt_k[0] = this->x(coord_ref_k);
        pt_k[1] = this->y(coord_ref_k);
        break;
      case 3:
        pt_k[0] = this->x(coord_ref_k);
        pt_k[1] = this->y(coord_ref_k);
        pt_k[2] = this->z(coord_ref_k);
        break;
    }

    residual = norm(pt_phy - pt_k);

    ++iter;
    if (iter > iter_max) {
      break;
    }
    if (abserr < 10. * std::numeric_limits< DataType >::epsilon())
      break;

  } // end newton

  LOG_DEBUG(1, "Inverse cell-trafo ended after "
                   << iter << " Newton iterations with residual = " << residual
                   << ", |x_k - x_{k-1}| = " << abserr);
  // Set values ...
  co_ref = ref_k;

  return residual < 10. * std::numeric_limits< DataType >::epsilon();
} // namespace doffem

//////////////////////////////////////////////////////////////////////Newton2Dto3D//////////////////////////////////////////////////////////////////////////////



template < class DataType, int DIM >
bool CellTransformation< DataType, DIM >::inverse_newton_2Dto3D( Vec<3, DataType> &co_phy, 
                                                                 Vec<2, DataType> &co_ref, 
                                                                 Vec<2, DataType> &co_ref_0) const 
{
/* needs to be debugged 
  // Initialisation

  Vec<3, DataType> pt_phy = co_phy;

  Vec<2, DataType> ref_k1, ref_k;
  ref_k = co_ref_0;

  const DataType tol_eps = 1.e3 * std::numeric_limits< DataType >::epsilon();

  // Some general parameters

  const int iter_max = 1000;
  int iter = 0;

  // Residual

  Vec<2, DataType> coord_ref_k = ref_k;
  Vec<3, DataType> pt_k;

   pt_k[0] = this->x(coord_ref_k);
   pt_k[1] = this->y(coord_ref_k);
   pt_k[2] = this->z(coord_ref_k);

  
  DataType residual = norm(pt_phy - pt_k);
  DataType abserr = norm(ref_k1 - ref_k);

  // Newton

  // Jacobian Matrix (grad G)
  Mat< 3, 2, DataType > G;
  // Pseudo-Inverse of the jacobian
  Mat< 2, 3, DataType > B;

  while (residual > 10. * std::numeric_limits< DataType >::epsilon()) 
  {
    G(0, 0) = this->x_x(coord_ref_k);
    G(0, 1) = this->x_y(coord_ref_k);

    G(1, 0) = this->y_x(coord_ref_k);
    G(1, 1) = this->y_y(coord_ref_k);

    G(2, 0) = this->z_x(coord_ref_k);
    G(2, 1) = this->z_y(coord_ref_k);

    //inverse of G^T * G
    
    Mat< 2, 2, DataType > res;

#ifndef NDEBUG
    const DataType detres = det(res);

    assert(det(G.transpose_me() * G) != 0.);
#endif

    inv(G.transpose_me() * G, res);

    B = res * G.transpose_me();

    // Armijo parameter

    const int iter_max_armijo = 500;
    int iter_armijo = 0;

    DataType residual_armijo = 2. * residual;

    DataType omega = 1.;

    const Vec<2, DataType> update_vec = B * (pt_phy - pt_k);

    // Start Armijo

    while ((iter_armijo <= iter_max_armijo) && (residual_armijo > residual)) {

      ref_k1 = ref_k + omega * update_vec;

      for (int d=0; d<2; ++d) 
      {
        if ((ref_k1[d] >= -tol_eps) && (ref_k1[d] <= tol_eps)) {
          ref_k1[d] = 0.;
        }
        if ((ref_k1[d] - 1. >= -tol_eps) && (ref_k1[d] - 1. <= tol_eps)) {
          ref_k1[d] = 1.;
        }
      }

      Vec<2, DataType> ref_check = ref_k1;

      while (!(this->contains_reference_point(ref_check))) {
        omega /= 2.0;
        ref_k1 = ref_k + omega * update_vec;

        for (int d=0; d<2; ++d) 
        {
          if ((ref_k1[d] >= -tol_eps) && (ref_k1[d] <= tol_eps)) {
            ref_k1[d] = 0.;
          }
          if ((ref_k1[d] - 1. >= -tol_eps) && (ref_k1[d] - 1. <= tol_eps)) {
            ref_k1[d] = 1.;
          }
        }
        ref_check = ref_k1;
      }

      Vec<2, DataType> coord_ref_k1 = ref_k1;

      Coord F_k1;

      F_k1[0] = this->x(coord_ref_k1);
      F_k1[1] = this->y(coord_ref_k1);
      F_k1[2] = this->z(coord_ref_k1);


      residual_armijo = norm(pt_phy - F_k1);

      ++iter_armijo;
      omega /= 2.;
    }

    abserr = norm(ref_k1 - ref_k);
    ref_k = ref_k1;

    coord_ref_k = ref_k;

    pt_k[0] = this->x(coord_ref_k);
    pt_k[1] = this->y(coord_ref_k);
    pt_k[2] = this->z(coord_ref_k);


    residual = norm(pt_phy - pt_k);

    ++iter;
    if (iter > iter_max) {
      break;
    }
    if (abserr < 10. * std::numeric_limits< DataType >::epsilon())
      break;

  } // end newton

  LOG_DEBUG(1, "Inverse cell-trafo ended after "
                   << iter << " Newton iterations with residual = " << residual
                   << ", |x_k - x_{k-1}| = " << abserr);
  // Set values ...
  co_ref = ref_k;

  return residual < 10. * std::numeric_limits< DataType >::epsilon();
*/
  return false;
} 

template < class DataType, int DIM >
bool CellTransformation< DataType, DIM >::contains_physical_point(
    const Coord &coord_phys, Coord& cr) const {

  bool found_ref_point = false;
  // Compute reference coordinates
  switch (DIM) {
  case 1: {
    found_ref_point = this->inverse(coord_phys, cr);
    LOG_DEBUG(2, "Physical Point " << coord_phys[0] << " has ref coords "
                                   << cr[0] << "? -> " << found_ref_point);
    break;
  }
  case 2: {
    found_ref_point = this->inverse(coord_phys, cr);
    LOG_DEBUG(2, "Physical Point (" << coord_phys[0] << ", " << coord_phys[1]
                                    << ") has ref coords " << cr[0] << ", "
                                    << cr[1] << ") ? -> " << found_ref_point);
    break;
  }
  case 3: {
    found_ref_point = this->inverse(coord_phys, cr);
    LOG_DEBUG(2, "Physical Point ("
                     << coord_phys[0] << ", " << coord_phys[1] << ", "
                     << coord_phys[2] << ") has ref coords " << cr[0] << ", "
                     << cr[1] << ", " << cr[2] << ") ? -> " << found_ref_point);
    break;
  }
  default: {
    throw "Invalid dimension!";
    break;
  }
  }

  if (!found_ref_point) {
    return false;
  }

  const bool contains_pt = this->contains_reference_point(cr);

  return contains_pt;
}

template<class DataType, int DIM>
void CellTransformation<DataType, DIM>::transform ( const Vec<DIM,DataType>& coord_ref, Vec<DIM, DataType>& coord_mapped ) const 
{
  switch (DIM)
  {
    case 1:
      coord_mapped[0] = this->x ( coord_ref );
      break;
    case 2:
      coord_mapped[0] = this->x ( coord_ref );
      coord_mapped[1] = this->y ( coord_ref );
      break;
    case 3:
      coord_mapped[0] = this->x ( coord_ref );
      coord_mapped[1] = this->y ( coord_ref );
      coord_mapped[2] = this->z ( coord_ref );
      break;
    default:
      assert(0);
      break;
  }
}

template<class DataType, int DIM>
void CellTransformation<DataType,DIM>::J ( const Vec<DIM,DataType>& coord, Mat<DIM, DIM, DataType>& mat) const 
{
  switch (DIM)
  {
    case 1:
      mat(0,0) = this->x_x(coord);
      break;
    case 2:
      mat(0,0) = this->x_x(coord);
      mat(0,1) = this->x_y(coord);
      mat(1,0) = this->y_x(coord);
      mat(1,1) = this->y_y(coord);
      break;
    case 3:
      mat(0,0) = this->x_x(coord);
      mat(0,1) = this->x_y(coord);
      mat(0,2) = this->x_z(coord);
      
      mat(1,0) = this->y_x(coord);
      mat(1,1) = this->y_y(coord);
      mat(1,2) = this->y_z(coord);

      mat(2,0) = this->z_x(coord);
      mat(2,1) = this->z_y(coord);
      mat(2,2) = this->z_z(coord);
      break;
    default:
      assert (0);
      break;
  }
}

template<class DataType, int DIM>
void CellTransformation<DataType,DIM>::J_and_detJ ( const Vec<DIM,DataType>& coord, Mat<DIM, DIM, DataType>& mat, DataType& detJ) const 
{
  switch (DIM)
  {
    case 1:
      mat(0,0) = this->x_x(coord);
      detJ = mat(0,0);
      break;
    case 2:
      mat(0,0) = this->x_x(coord);
      mat(0,1) = this->x_y(coord);
      mat(1,0) = this->y_x(coord);
      mat(1,1) = this->y_y(coord);
      detJ = mat(0,0) * mat(1,1) - mat(0,1) * mat(1,0);
      break;
    case 3:
      mat(0,0) = this->x_x(coord);
      mat(0,1) = this->x_y(coord);
      mat(0,2) = this->x_z(coord);
      
      mat(1,0) = this->y_x(coord);
      mat(1,1) = this->y_y(coord);
      mat(1,2) = this->y_z(coord);

      mat(2,0) = this->z_x(coord);
      mat(2,1) = this->z_y(coord);
      mat(2,2) = this->z_z(coord);
      
      detJ = mat(0,0) * (mat(1,1) * mat(2,2) - mat(1,2) *  mat(2,1))
           - mat(0,1) * (mat(1,0) * mat(2,2) - mat(1,2) *  mat(2,0))
           + mat(0,2) * (mat(1,0) * mat(2,1) - mat(1,1) *  mat(2,0));
      break;
    default:
      assert (0);
      break;
  }
}

template<class DataType, int DIM>
void CellTransformation<DataType, DIM>::H ( const Vec<DIM,DataType>& coord, size_t d, Mat<DIM, DIM, DataType>& mat) const 
{
  switch (DIM)
  {
    case 1:
      switch (d)
      {
        case 0:
          mat(0,0) = this->x_xx(coord);
          break;
        default:
          assert(0);
          break;
      }
      break;
    case 2:
      switch (d)
      {
        case 0:
          mat(0,0) = this->x_xx(coord);
          mat(0,1) = this->x_xy(coord);
          mat(1,0) = mat(0,1);
          mat(1,1) = this->x_yy(coord);
          break;
        case 1:
          mat(0,0) = this->y_xx(coord);
          mat(0,1) = this->y_xy(coord);
          mat(1,0) = mat(0,1);
          mat(1,1) = this->y_yy(coord);
          break;
        default:
          assert(0);
          break;
      }
      break;
    case 3:
      switch (d)
      {
        case 0:
          mat(0,0) = this->x_xx(coord);
          mat(0,1) = this->x_xy(coord);
          mat(0,2) = this->x_xz(coord);
          mat(1,0) = mat(0,1);
          mat(1,1) = this->x_yy(coord);
          mat(1,2) = this->x_yz(coord);
          mat(2,0) = mat(0,2);
          mat(2,1) = mat(1,2);
          mat(2,2) = this->x_zz(coord);
          break;
        case 1:
          mat(0,0) = this->y_xx(coord);
          mat(0,1) = this->y_xy(coord);
          mat(0,2) = this->y_xz(coord);
          mat(1,0) = mat(0,1);
          mat(1,1) = this->y_yy(coord);
          mat(1,2) = this->y_yz(coord);
          mat(2,0) = mat(0,2);
          mat(2,1) = mat(1,2);
          mat(2,2) = this->y_zz(coord);
          break;
        case 2:
          mat(0,0) = this->z_xx(coord);
          mat(0,1) = this->z_xy(coord);
          mat(0,2) = this->z_xz(coord);
          mat(1,0) = mat(0,1);
          mat(1,1) = this->z_yy(coord);
          mat(1,2) = this->z_yz(coord);
          mat(2,0) = mat(0,2);
          mat(2,1) = mat(1,2);
          mat(2,2) = this->z_zz(coord);
          break;
        default:
          assert(0);
          break;
      }
      break;
    default:
      assert (0);
      break;
  }
}

template<class DataType, int DIM>
DataType CellTransformation<DataType, DIM>::detJ ( const Vec<DIM,DataType>& pt) const 
{
  switch (DIM)
  {
    case 1:
      return this->x_x(pt);
    case 2:
      return this->x_x(pt) * this->y_y(pt) - this->x_y(pt) * this->y_x(pt);
    case 3:
      return this->x_x(pt) * (this->y_y(pt) * this->z_z(pt) - this->y_z(pt) * this->z_y(pt)) 
           - this->x_y(pt) * (this->y_x(pt) * this->z_z(pt) - this->y_z(pt) * this->z_x(pt))
           + this->x_z(pt) * (this->y_x(pt) * this->z_y(pt) - this->y_y(pt) * this->z_x(pt));
    default:
      assert(0);
      break;
  }
}

template<class DataType, int DIM>
void CellTransformation<DataType,DIM>::grad_detJ ( const Vec<DIM,DataType>& coord, Vec<DIM, DataType>& grad) const 
{
  switch (DIM)
  {
    case 1:
      grad[0] = this->detJ_x(coord);
      break;
    case 2:
      grad[0] = this->detJ_x(coord);
      grad[1] = this->detJ_y(coord);
      break;
    case 3:
      grad[0] = this->detJ_x(coord);
      grad[1] = this->detJ_y(coord);
      grad[2] = this->detJ_z(coord);
      break;
    default:
      assert (0);
      break;
  }
}

template<class DataType, int DIM>
void CellTransformation<DataType,DIM>::grad_inv_detJ ( const Vec<DIM,DataType>& coord, Vec<DIM, DataType>& grad) const 
{
  const DataType detJ = this->detJ(coord);
  const DataType inv_detJ_2 = -1. / (detJ * detJ);
  
  switch (DIM)
  {
    case 1:
      grad[0] = inv_detJ_2 * this->detJ_x(coord);
      break;
    case 2:
      grad[0] = inv_detJ_2 * this->detJ_x(coord);
      grad[1] = inv_detJ_2 * this->detJ_y(coord);
      break;
    case 3:
      grad[0] = inv_detJ_2 * this->detJ_x(coord);
      grad[1] = inv_detJ_2 * this->detJ_y(coord);
      grad[2] = inv_detJ_2 * this->detJ_z(coord);
      break;
    default:
      assert (0);
      break;
  }
}

template<class DataType, int DIM>
void CellTransformation<DataType,DIM>::hessian_detJ ( const Vec<DIM,DataType>& coord, Mat<DIM, DIM, DataType>& mat) const 
{
  switch (DIM)
  {
    case 1:
      mat(0,0) = this->detJ_xx(coord);
      break;
    case 2:
      mat(0,0) = this->detJ_xx(coord);
      mat(0,1) = this->detJ_xy(coord);
      mat(1,0) = mat(0,1);
      mat(1,1) = this->detJ_yy(coord);
      break;
    case 3:
      mat(0,0) = this->detJ_xx(coord);
      mat(0,1) = this->detJ_xy(coord);
      mat(0,2) = this->detJ_xz(coord);
      
      mat(1,0) = mat(0,1);
      mat(1,1) = this->detJ_yy(coord);
      mat(1,2) = this->detJ_yz(coord);

      mat(2,0) = mat(0,2);
      mat(2,1) = mat(1,2);
      mat(2,2) = this->detJ_zz(coord);
      break;
    default:
      assert (0);
      break;
  }
}

template<class DataType, int DIM>
DataType CellTransformation<DataType, DIM>::detJ_x ( const Vec<DIM,DataType>& pt) const 
{
  if (this->order_ <= 1)
  {
    return 0.;
  }
  switch (DIM)
  {
    case 1:
      return this->x_xx(pt);
    case 2:
      return this->x_xx(pt) * this->y_y(pt) + this->x_x(pt) * this->y_xy(pt)  
           - this->x_xy(pt) * this->y_x(pt) - this->x_y(pt) * this->y_xx(pt); 
    case 3:
      return  this->x_xx(pt) * (this->y_y(pt)  * this->z_z(pt) - this->y_z(pt) * this->z_y(pt))
            + this->x_x(pt)  * (this->y_xy(pt) * this->z_z(pt) + this->y_y(pt) * this->z_xz(pt) 
                              - this->y_xz(pt) * this->z_y(pt) - this->y_z(pt) * this->z_xy(pt))
            - this->x_xy(pt) * (this->y_x(pt)  * this->z_z(pt) - this->y_z(pt) * this->z_x(pt))
            - this->x_y(pt)  * (this->y_xx(pt) * this->z_z(pt) + this->y_x(pt) * this->z_xz(pt) 
                              - this->y_xz(pt) * this->z_x(pt) - this->y_z(pt) * this->z_xx(pt))
            + this->x_xz(pt) * (this->y_x(pt)  * this->z_y(pt) - this->y_y(pt) * this->z_x(pt))
            + this->x_z(pt)  * (this->y_xx(pt) * this->z_y(pt) + this->y_x(pt) * this->z_xy(pt)
                              - this->y_xy(pt) * this->z_x(pt) - this->y_y(pt) * this->z_xx(pt));
    default:
      assert(0);
      break;
  }
}

template<class DataType, int DIM>
DataType CellTransformation<DataType, DIM>::detJ_y ( const Vec<DIM,DataType>& pt) const 
{
  if (this->order_ <= 1)
  {
    return 0.;
  }
  switch (DIM)
  {
    case 1:
      return 0.;
    case 2:
      return this->x_xy(pt) * this->y_y(pt) + this->x_x(pt) * this->y_yy(pt)  
           - this->x_yy(pt) * this->y_x(pt) - this->x_y(pt) * this->y_xy(pt); 
    case 3:
      return  this->x_xy(pt) * (this->y_y(pt)  * this->z_z(pt) - this->y_z(pt) * this->z_y(pt))
            + this->x_x(pt)  * (this->y_yy(pt) * this->z_z(pt) + this->y_y(pt) * this->z_yz(pt) 
                              - this->y_yz(pt) * this->z_y(pt) - this->y_z(pt) * this->z_yy(pt))
            - this->x_yy(pt) * (this->y_x(pt)  * this->z_z(pt) - this->y_z(pt) * this->z_x(pt))
            - this->x_y(pt)  * (this->y_xy(pt) * this->z_z(pt) + this->y_x(pt) * this->z_yz(pt) 
                              - this->y_yz(pt) * this->z_x(pt) - this->y_z(pt) * this->z_xy(pt))
            + this->x_yz(pt) * (this->y_x(pt)  * this->z_y(pt) - this->y_y(pt) * this->z_x(pt))
            + this->x_z(pt)  * (this->y_xy(pt) * this->z_y(pt) + this->y_x(pt) * this->z_yy(pt)
                              - this->y_yy(pt) * this->z_x(pt) - this->y_y(pt) * this->z_xy(pt));  
    default:
      assert(0);
      break;
  }
}

template<class DataType, int DIM>
DataType CellTransformation<DataType, DIM>::detJ_z ( const Vec<DIM,DataType>& pt) const 
{
  if (this->order_ <= 1)
  {
    return 0.;
  }
  switch (DIM)
  {
    case 1:
      return 0.;
    case 2:
      return 0.;
    case 3:
      return  this->x_xz(pt) * (this->y_y(pt)  * this->z_z(pt) - this->y_z(pt) * this->z_y(pt))
            + this->x_x(pt)  * (this->y_yz(pt) * this->z_z(pt) + this->y_y(pt) * this->z_zz(pt) 
                              - this->y_zz(pt) * this->z_y(pt) - this->y_z(pt) * this->z_yz(pt))
            - this->x_yz(pt) * (this->y_x(pt)  * this->z_z(pt) - this->y_z(pt) * this->z_x(pt))
            - this->x_y(pt)  * (this->y_xz(pt) * this->z_z(pt) + this->y_x(pt) * this->z_zz(pt) 
                              - this->y_zz(pt) * this->z_x(pt) - this->y_z(pt) * this->z_xz(pt))
            + this->x_zz(pt) * (this->y_x(pt)  * this->z_y(pt) - this->y_y(pt) * this->z_x(pt))
            + this->x_z(pt)  * (this->y_xz(pt) * this->z_y(pt) + this->y_x(pt) * this->z_yz(pt)
                              - this->y_yz(pt) * this->z_x(pt) - this->y_y(pt) * this->z_xz(pt));                                            
    default:
      assert(0);
      break;
  }
}

template<class DataType, int DIM>
DataType CellTransformation<DataType, DIM>::detJ_xx ( const Vec<DIM,DataType>& pt) const 
{
  if (this->order_ <= 2)
  {
    return 0.;
  }
  switch (DIM)
  {
    case 1:
      return  this->x_xxx(pt);
    case 2:
      return  this->x_xxx(pt) * this->y_y(pt)  + this->x_xx(pt) * this->y_xy(pt) 
            + this->x_xx(pt)  * this->y_xy(pt) + this->x_x(pt)  * this->y_xxy(pt)  
            - this->x_xxy(pt) * this->y_x(pt)  - this->x_xy(pt) * this->y_xx(pt) 
            - this->x_xy(pt)  * this->y_xx(pt) - this->x_y(pt)  * this->y_xxx(pt) ;  
    case 3:
      return       this->x_xxx(pt) * (this->y_y(pt)   * this->z_z(pt)  - this->y_z(pt)  * this->z_y(pt))
            + 2. * this->x_xx(pt)  * (this->y_xy(pt)  * this->z_z(pt)  + this->y_y(pt)  * this->z_xz(pt) 
                                    - this->y_xz(pt)  * this->z_y(pt)  - this->y_z(pt)  * this->z_xy(pt))
            +      this->x_x(pt)   * (this->y_xxy(pt) * this->z_z(pt)  + this->y_xy(pt) * this->z_xz(pt) 
                                    + this->y_xy(pt)  * this->z_xz(pt) + this->y_y(pt)  * this->z_xxz(pt)
                                    - this->y_xxz(pt) * this->z_y(pt)  - this->y_xz(pt) * this->z_xy(pt) 
                                    - this->y_xz(pt)  * this->z_xy(pt) - this->y_z(pt)  * this->z_xxy(pt))
            -      this->x_xxy(pt) * (this->y_x(pt)   * this->z_z(pt)  - this->y_z(pt)  * this->z_x(pt))
            - 2. * this->x_xy(pt)  * (this->y_xx(pt)  * this->z_z(pt)  + this->y_x(pt)  * this->z_xz(pt) 
                                    - this->y_xz(pt)  * this->z_x(pt)  - this->y_z(pt)  * this->z_xx(pt))
            -      this->x_y(pt)   * (this->y_xxx(pt) * this->z_z(pt)  + this->y_xx(pt) * this->z_xz(pt) 
                                    + this->y_xx(pt)  * this->z_xz(pt) + this->y_x(pt)  * this->z_xxz(pt)  
                                    - this->y_xxz(pt) * this->z_x(pt)  - this->y_xz(pt) * this->z_xx(pt) 
                                    - this->y_xz(pt)  * this->z_xx(pt) - this->y_z(pt)  * this->z_xxx(pt))
            +      this->x_xxz(pt) * (this->y_x(pt)   * this->z_y(pt)  - this->y_y(pt)  * this->z_x(pt))
            + 2. * this->x_xz(pt)  * (this->y_xx(pt)  * this->z_y(pt)  + this->y_x(pt)  * this->z_xy(pt)
                                    - this->y_xy(pt)  * this->z_x(pt)  - this->y_y(pt)  * this->z_xx(pt))
            +      this->x_z(pt)   * (this->y_xxx(pt) * this->z_y(pt)  + this->y_xx(pt) * this->z_xy(pt)
                                    + this->y_xx(pt)  * this->z_xy(pt) + this->y_x(pt)  * this->z_xxy(pt)
                                    - this->y_xxy(pt) * this->z_x(pt)  - this->y_xy(pt) * this->z_xx(pt)
                                    - this->y_xy(pt)  * this->z_xx(pt) - this->y_y(pt)  * this->z_xxx(pt));
    default:
      assert(0);
      break;
  }
}

template<class DataType, int DIM>
DataType CellTransformation<DataType, DIM>::detJ_xy ( const Vec<DIM,DataType>& pt) const 
{
  if (this->order_ <= 2)
  {
    return 0.;
  }
  switch (DIM)
  {
    case 1:
      return 0.;
    case 2:
      return  this->x_xxy(pt) * this->y_y(pt)  + this->x_xx(pt) * this->y_yy(pt) 
          + this->x_xy(pt)  * this->y_xy(pt) + this->x_x(pt)  * this->y_xyy(pt)  
          - this->x_xyy(pt) * this->y_x(pt)  - this->x_xy(pt) * this->y_xy(pt) 
          - this->x_yy(pt)  * this->y_xx(pt) - this->x_y(pt)  * this->y_xxy(pt);   
    case 3:
      return  this->x_xxy(pt) * (this->y_y(pt)   * this->z_z(pt)  - this->y_z(pt)  * this->z_y(pt))
            + this->x_xy(pt)  * (this->y_yy(pt)  * this->z_z(pt)  + this->y_y(pt)  * this->z_yz(pt)
                               - this->y_yz(pt)  * this->z_y(pt)  - this->y_z(pt)  * this->z_yy(pt))
            + this->x_xy(pt)  * (this->y_xy(pt)  * this->z_z(pt)  + this->y_y(pt)  * this->z_xz(pt) 
                               - this->y_xz(pt)  * this->z_y(pt)  - this->y_z(pt)  * this->z_xy(pt))
            + this->x_x(pt)   * (this->y_xyy(pt) * this->z_z(pt)  + this->y_xy(pt) * this->z_yz(pt) 
                               + this->y_yy(pt)  * this->z_xz(pt) + this->y_y(pt)  * this->z_xyz(pt) 
                               - this->y_xyz(pt) * this->z_y(pt)  - this->y_xz(pt) * this->z_yy(pt) 
                               - this->y_yz(pt)  * this->z_xy(pt) - this->y_z(pt)  * this->z_xyy(pt))
            - this->x_xyy(pt) * (this->y_x(pt)   * this->z_z(pt)  - this->y_z(pt)  * this->z_x(pt))
            - this->x_xy(pt)  * (this->y_xy(pt)  * this->z_z(pt)  + this->y_x(pt)  * this->z_yz(pt) 
                               - this->y_yz(pt)  * this->z_x(pt)  - this->y_z(pt)  * this->z_xy(pt))
            - this->x_yy(pt)  * (this->y_xx(pt)  * this->z_z(pt)  + this->y_x(pt)  * this->z_xz(pt) 
                               - this->y_xz(pt)  * this->z_x(pt)  - this->y_z(pt)  * this->z_xx(pt))        
            - this->x_y(pt)   * (this->y_xxy(pt) * this->z_z(pt)  + this->y_xx(pt) * this->z_yz(pt) 
                               + this->y_xy(pt)  * this->z_xz(pt) + this->y_x(pt)  * this->z_xyz(pt) 
                               - this->y_xyz(pt) * this->z_x(pt)  - this->y_xz(pt) * this->z_xy(pt) 
                               - this->y_yz(pt)  * this->z_xx(pt) - this->y_z(pt)  * this->z_xxy(pt))
            + this->x_xyz(pt) * (this->y_x(pt)   * this->z_y(pt)  - this->y_y(pt)  * this->z_x(pt))
            + this->x_xz(pt)  * (this->y_xy(pt)  * this->z_y(pt)  + this->y_x(pt)  * this->z_yy(pt) 
                               - this->y_yy(pt)  * this->z_x(pt)  - this->y_y(pt)  * this->z_xy(pt))
            + this->x_yz(pt)  * (this->y_xx(pt)  * this->z_y(pt)  + this->y_x(pt)  * this->z_xy(pt)
                               - this->y_xy(pt)  * this->z_x(pt)  - this->y_y(pt)  * this->z_xx(pt))
            + this->x_z(pt)   * (this->y_xxy(pt) * this->z_y(pt)  + this->y_xx(pt) * this->z_yy(pt) 
                               + this->y_xy(pt)  * this->z_xy(pt) + this->y_x(pt)  * this->z_xyy(pt)
                               - this->y_xyy(pt) * this->z_x(pt)  - this->y_xy(pt) * this->z_xy(pt) 
                               - this->y_yy(pt)  * this->z_xx(pt) - this->y_y(pt)  * this->z_xxy(pt)); 
    default:
      assert(0);
      break;
  }
}

template<class DataType, int DIM>
DataType CellTransformation<DataType, DIM>::detJ_xz ( const Vec<DIM,DataType>& pt) const 
{
  if (this->order_ <= 2)
  {
    return 0.;
  }
  switch (DIM)
  {
    case 1:
      return 0.;
    case 2:
      return 0.;
    case 3:
      return  this->x_xxz(pt) * (this->y_y(pt)   * this->z_z(pt)  - this->y_z(pt)  * this->z_y(pt))
            + this->x_xx(pt)  * (this->y_yz(pt)  * this->z_z(pt)  + this->y_y(pt)  * this->z_zz(pt) 
                               - this->y_zz(pt)  * this->z_y(pt)  - this->y_z(pt)  * this->z_yz(pt))
            + this->x_xz(pt)  * (this->y_xy(pt)  * this->z_z(pt)  + this->y_y(pt)  * this->z_xz(pt) 
                               - this->y_xz(pt)  * this->z_y(pt)  - this->y_z(pt)  * this->z_xy(pt))
            + this->x_x(pt)   * (this->y_xyz(pt) * this->z_z(pt)  + this->y_xy(pt) * this->z_zz(pt) 
                               + this->y_yz(pt)  * this->z_xz(pt) + this->y_y(pt)  * this->z_xzz(pt) 
                               - this->y_xzz(pt) * this->z_y(pt)  - this->y_xz(pt) * this->z_yz(pt)
                               - this->y_zz(pt)  * this->z_xy(pt) - this->y_z(pt)  * this->z_xyz(pt))
            - this->x_xyz(pt) * (this->y_x(pt)   * this->z_z(pt)  - this->y_z(pt)  * this->z_x(pt))
            - this->x_xy(pt)  * (this->y_xz(pt)  * this->z_z(pt)  + this->y_x(pt)  * this->z_zz(pt) 
                               - this->y_zz(pt)  * this->z_x(pt)  - this->y_z(pt)  * this->z_xz(pt))
            - this->x_yz(pt)  * (this->y_xx(pt)  * this->z_z(pt)  + this->y_x(pt)  * this->z_xz(pt) 
                               - this->y_xz(pt)  * this->z_x(pt)  - this->y_z(pt)  * this->z_xx(pt))
            - this->x_y(pt)   * (this->y_xxz(pt) * this->z_z(pt)  + this->y_xx(pt) * this->z_zz(pt) 
                               + this->y_xz(pt)  * this->z_xz(pt) + this->y_x(pt)  * this->z_xzz(pt) 
                               - this->y_xzz(pt) * this->z_x(pt)  - this->y_xz(pt) * this->z_xz(pt) 
                               - this->y_zz(pt)  * this->z_xx(pt) - this->y_z(pt)  * this->z_xxz(pt))
            + this->x_xzz(pt) * (this->y_x(pt)   * this->z_y(pt)  - this->y_y(pt)  * this->z_x(pt))
            + this->x_xz(pt)  * (this->y_xz(pt)  * this->z_y(pt)  + this->y_x(pt)  * this->z_yz(pt) 
                               - this->y_yz(pt)  * this->z_x(pt)  - this->y_y(pt)  * this->z_xz(pt))
            + this->x_zz(pt)  * (this->y_xx(pt)  * this->z_y(pt)  + this->y_x(pt)  * this->z_xy(pt)
                               - this->y_xy(pt)  * this->z_x(pt)  - this->y_y(pt)  * this->z_xx(pt))                                            
            + this->x_z(pt)   * (this->y_xxz(pt) * this->z_y(pt)  + this->y_xx(pt) * this->z_yz(pt) 
                               + this->y_xz(pt)  * this->z_xy(pt) + this->y_x(pt)  * this->z_xyz(pt)
                               - this->y_xyz(pt) * this->z_x(pt)  - this->y_xy(pt) * this->z_xz(pt) 
                               - this->y_yz(pt)  * this->z_xx(pt) - this->y_y(pt)  * this->z_xxz(pt));  
    default:
      assert(0);
      break;
  }
}

template<class DataType, int DIM>
DataType CellTransformation<DataType, DIM>::detJ_yy ( const Vec<DIM,DataType>& pt) const 
{
  if (this->order_ <= 2)
  {
    return 0.;
  }
  switch (DIM)
  {
    case 1:
      return 0.;
    case 2:
      return this->x_xyy(pt) * this->y_y(pt)  + this->x_xy(pt) * this->y_yy(pt) 
           + this->x_xy(pt)  * this->y_yy(pt) + this->x_x(pt)  * this->y_yyy(pt)  
           - this->x_yyy(pt) * this->y_x(pt)  - this->x_yy(pt) * this->y_xy(pt) 
           - this->x_yy(pt)  * this->y_xy(pt) - this->x_y(pt)  * this->y_xyy(pt); 
    case 3:
      return  this->x_xyy(pt) * (this->y_y(pt)   * this->z_z(pt)  - this->y_z(pt)  * this->z_y(pt))
            + this->x_xy(pt)  * (this->y_yy(pt)  * this->z_z(pt)  + this->y_y(pt)  * this->z_yz(pt) 
                               - this->y_yz(pt)  * this->z_y(pt)  - this->y_z(pt)  * this->z_yy(pt))
            + this->x_xy(pt)  * (this->y_yy(pt)  * this->z_z(pt)  + this->y_y(pt)  * this->z_yz(pt) 
                               - this->y_yz(pt)  * this->z_y(pt)  - this->y_z(pt)  * this->z_yy(pt))
            + this->x_x(pt)   * (this->y_yyy(pt) * this->z_z(pt)  + this->y_yy(pt) * this->z_yz(pt) 
                               + this->y_yy(pt)  * this->z_yz(pt) + this->y_y(pt)  * this->z_yyz(pt) 
                               - this->y_yyz(pt) * this->z_y(pt)  - this->y_yz(pt) * this->z_yy(pt) 
                               - this->y_yz(pt)  * this->z_yy(pt) - this->y_z(pt)  * this->z_yyy(pt))
            - this->x_yyy(pt) * (this->y_x(pt)   * this->z_z(pt)  - this->y_z(pt)  * this->z_x(pt))
            - this->x_yy(pt)  * (this->y_xy(pt)  * this->z_z(pt)  + this->y_x(pt)  * this->z_yz(pt) 
                               - this->y_yz(pt)  * this->z_x(pt)  - this->y_z(pt)  * this->z_xy(pt))
            - this->x_yy(pt)  * (this->y_xy(pt)  * this->z_z(pt)  + this->y_x(pt)  * this->z_yz(pt) 
                               - this->y_yz(pt)  * this->z_x(pt)  - this->y_z(pt)  * this->z_xy(pt))
            - this->x_y(pt)   * (this->y_xyy(pt) * this->z_z(pt)  + this->y_xy(pt) * this->z_yz(pt) 
                               + this->y_xy(pt)  * this->z_yz(pt) + this->y_x(pt)  * this->z_yyz(pt) 
                               - this->y_yyz(pt) * this->z_x(pt)  - this->y_yz(pt) * this->z_xy(pt) 
                               - this->y_yz(pt)  * this->z_xy(pt) - this->y_z(pt)  * this->z_xyy(pt))
            + this->x_yyz(pt) * (this->y_x(pt)   * this->z_y(pt)  - this->y_y(pt)  * this->z_x(pt))
            + this->x_yz(pt)  * (this->y_xy(pt)  * this->z_y(pt)  + this->y_x(pt)  * this->z_yy(pt) 
                               - this->y_yy(pt)  * this->z_x(pt)  - this->y_y(pt)  * this->z_xy(pt))
            + this->x_yz(pt)  * (this->y_xy(pt)  * this->z_y(pt)  + this->y_x(pt)  * this->z_yy(pt)
                               - this->y_yy(pt)  * this->z_x(pt)  - this->y_y(pt)  * this->z_xy(pt)) 
            + this->x_z(pt)   * (this->y_xyy(pt) * this->z_y(pt)  + this->y_xy(pt) * this->z_yy(pt) 
                               + this->y_xy(pt)  * this->z_yy(pt) + this->y_x(pt)  * this->z_yyy(pt)
                               - this->y_yyy(pt) * this->z_x(pt)  - this->y_yy(pt) * this->z_xy(pt) 
                               - this->y_yy(pt)  * this->z_xy(pt) - this->y_y(pt)  * this->z_xyy(pt));
    default:
      assert(0);
      break;
  }
}

template<class DataType, int DIM>
DataType CellTransformation<DataType, DIM>::detJ_yz ( const Vec<DIM,DataType>& pt) const 
{
  if (this->order_ <= 2)
  {
    return 0.;
  }
  switch (DIM)
  {
    case 1:
      return 0.;
    case 2:
      return 0.;
    case 3:
      return  this->x_xyz(pt) * (this->y_y(pt)   * this->z_z(pt)  - this->y_z(pt)  * this->z_y(pt))
            + this->x_xy(pt)  * (this->y_yz(pt)  * this->z_z(pt)  + this->y_y(pt)  * this->z_zz(pt) 
                               - this->y_zz(pt)  * this->z_y(pt)  - this->y_z(pt)  * this->z_yz(pt))
            + this->x_xz(pt)  * (this->y_yy(pt)  * this->z_z(pt)  + this->y_y(pt)  * this->z_yz(pt) 
                               - this->y_yz(pt)  * this->z_y(pt)  - this->y_z(pt)  * this->z_yy(pt))
            + this->x_x(pt)   * (this->y_yyz(pt) * this->z_z(pt)  + this->y_yy(pt) * this->z_zz(pt) 
                               + this->y_yz(pt)  * this->z_yz(pt) + this->y_y(pt)  * this->z_yzz(pt) 
                               - this->y_yzz(pt) * this->z_y(pt)  - this->y_yz(pt) * this->z_yz(pt) 
                               - this->y_zz(pt)  * this->z_yy(pt) - this->y_z(pt)  * this->z_yyz(pt))
            - this->x_yyz(pt) * (this->y_x(pt)   * this->z_z(pt)  - this->y_z(pt)  * this->z_x(pt))
            - this->x_yy(pt)  * (this->y_xz(pt)  * this->z_z(pt)  + this->y_x(pt)  * this->z_zz(pt) 
                               - this->y_zz(pt)  * this->z_x(pt)  - this->y_z(pt)  * this->z_xz(pt))
            - this->x_yz(pt)  * (this->y_xy(pt)  * this->z_z(pt)  + this->y_x(pt)  * this->z_yz(pt) 
                               - this->y_yz(pt)  * this->z_x(pt)  - this->y_z(pt)  * this->z_xy(pt))
            - this->x_y(pt)   * (this->y_xyz(pt) * this->z_z(pt)  + this->y_xy(pt) * this->z_zz(pt) 
                               + this->y_xz(pt)  * this->z_yz(pt) + this->y_x(pt)  * this->z_yzz(pt) 
                               - this->y_yzz(pt) * this->z_x(pt)  - this->y_yz(pt) * this->z_xz(pt) 
                               - this->y_zz(pt)  * this->z_xy(pt) - this->y_z(pt)  * this->z_xyz(pt))
            + this->x_yzz(pt) * (this->y_x(pt)   * this->z_y(pt)  - this->y_y(pt)  * this->z_x(pt))
            + this->x_yz(pt)  * (this->y_xz(pt)  * this->z_y(pt)  + this->y_x(pt)  * this->z_yz(pt) 
                               - this->y_yz(pt)  * this->z_x(pt)  - this->y_y(pt)  * this->z_xz(pt))
            + this->x_zz(pt)  * (this->y_xy(pt)  * this->z_y(pt)  + this->y_x(pt)  * this->z_yy(pt)
                               - this->y_yy(pt)  * this->z_x(pt)  - this->y_y(pt)  * this->z_xy(pt))                                            
            + this->x_z(pt)   * (this->y_xyz(pt) * this->z_y(pt)  + this->y_xy(pt) * this->z_yz(pt) 
                               + this->y_xz(pt)  * this->z_yy(pt) + this->y_x(pt)  * this->z_yyz(pt)
                               - this->y_yyz(pt) * this->z_x(pt)  - this->y_yy(pt) * this->z_xz(pt) 
                               - this->y_yz(pt)  * this->z_xy(pt) - this->y_y(pt)  * this->z_xyz(pt));
    default:
      assert(0);
      break;
  }
}

template<class DataType, int DIM>
DataType CellTransformation<DataType, DIM>::detJ_zz ( const Vec<DIM,DataType>& pt) const 
{
  if (this->order_ <= 2)
  {
    return 0.;
  }
  switch (DIM)
  {
    case 1:
      return 0.;
    case 2:
      return 0.;
    case 3:
      return  this->x_xzz(pt) * (this->y_y(pt)   * this->z_z(pt)  - this->y_z(pt)  * this->z_y(pt))
            +  this->x_xz(pt) * (this->y_yz(pt)  * this->z_z(pt)  + this->y_y(pt)  * this->z_zz(pt) 
                               - this->y_zz(pt)  * this->z_y(pt)  - this->y_z(pt)  * this->z_yz(pt))
            + this->x_xz(pt)  * (this->y_yz(pt)  * this->z_z(pt)  + this->y_y(pt)  * this->z_zz(pt) 
                               - this->y_zz(pt)  * this->z_y(pt)  - this->y_z(pt)  * this->z_yz(pt))
            + this->x_x(pt)   * (this->y_yzz(pt) * this->z_z(pt)  + this->y_yz(pt) * this->z_zz(pt) 
                               + this->y_yz(pt)  * this->z_zz(pt) + this->y_y(pt)  * this->z_zzz(pt) 
                               - this->y_zzz(pt) * this->z_y(pt)  - this->y_zz(pt) * this->z_yz(pt) 
                               - this->y_zz(pt)  * this->z_yz(pt) - this->y_z(pt)  * this->z_yzz(pt))
            - this->x_yzz(pt) * (this->y_x(pt)   * this->z_z(pt)  - this->y_z(pt)  * this->z_x(pt))
            - this->x_yz(pt)  * (this->y_xz(pt)  * this->z_z(pt)  + this->y_x(pt)  * this->z_zz(pt) 
                               - this->y_zz(pt)  * this->z_x(pt)  - this->y_z(pt)  * this->z_xz(pt))
            - this->x_yz(pt)  * (this->y_xz(pt)  * this->z_z(pt)  + this->y_x(pt)  * this->z_zz(pt) 
                               - this->y_zz(pt)  * this->z_x(pt)  - this->y_z(pt)  * this->z_xz(pt))
            - this->x_y(pt)   * (this->y_xzz(pt) * this->z_z(pt)  + this->y_xz(pt) * this->z_zz(pt)
                               + this->y_xz(pt)  * this->z_zz(pt) + this->y_x(pt)  * this->z_zzz(pt) 
                               - this->y_zzz(pt) * this->z_x(pt)  - this->y_zz(pt) * this->z_xz(pt) 
                               - this->y_zz(pt)  * this->z_xz(pt) - this->y_z(pt)  * this->z_xzz(pt))
            + this->x_zzz(pt) * (this->y_x(pt)   * this->z_y(pt)  - this->y_y(pt)  * this->z_x(pt))
            + this->x_zz(pt)  * (this->y_xz(pt)  * this->z_y(pt)  + this->y_x(pt)  * this->z_yz(pt) 
                               - this->y_yz(pt)  * this->z_x(pt)  - this->y_y(pt)  * this->z_xz(pt))
            + this->x_zz(pt)  * (this->y_xz(pt)  * this->z_y(pt)  + this->y_x(pt)  * this->z_yz(pt)
                               - this->y_yz(pt)  * this->z_x(pt)  - this->y_y(pt)  * this->z_xz(pt))                                          
            + this->x_z(pt)   * (this->y_xzz(pt) * this->z_y(pt)  + this->y_xz(pt) * this->z_yz(pt) 
                               + this->y_xz(pt)  * this->z_yz(pt) + this->y_x(pt)  * this->z_yzz(pt)
                               - this->y_yzz(pt) * this->z_x(pt)  - this->y_yz(pt) * this->z_xz(pt) 
                               - this->y_yz(pt)  * this->z_xz(pt) - this->y_y(pt)  * this->z_xzz(pt));
    default:
      assert(0);
      break;
  }
}

} // namespace doffem
} // namespace hiflow

#endif
