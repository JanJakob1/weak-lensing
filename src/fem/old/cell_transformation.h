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

#ifndef __FEM_CELL_TRANSFORMATION_H_
#define __FEM_CELL_TRANSFORMATION_H_

#include <cassert>
#include <vector>

#include "common/macros.h"
#include "dof/dof_fem_types.h"
#include <limits>

namespace hiflow {
namespace doffem {

///
/// \class CellTransformation cell_transformation.h
/// \brief Ancestor class of all transformation mappings from reference to
/// physical cells \author Michael Schick<br>Martin Baumann
///

template < class DataType > class CellTransformation {
public:
  typedef std::vector< DataType > Coord;

  /// Use this constructor which needs the geometrical dimension as input
  explicit CellTransformation(int gdim);

  virtual ~CellTransformation() {}

  /// Reinitialization of the transformation via coordinates of physical cell
  virtual void reinit(const Coord &coord_vtx);

  Coord get_coordinates() const { return this->coord_vtx_; }

  /// return vertex coordinates of reference cell
  virtual std::vector< Coord > get_reference_coordinates() const = 0;

  /// \brief Given physical cell coordinates in 1D,
  ///        this routine computes the corresponding reference cell coordinates
  /// @return true, if inverse computation was successful
  virtual bool inverse(DataType x_phy, DataType &x_ref) const = 0;

  /// \brief Given physical cell coordinates in 2D,
  ///        this routine computes the corresponding reference cell coordinates
  /// @return true, if inverse computation was successful
  virtual bool inverse(DataType x_phy, DataType y_phy, DataType &x_ref,
                       DataType &y_ref) const = 0;

  /// \brief Given physical cell coordinates in 3D,
  ///        this routine computes the corresponding reference cell coordinates
  /// @return true, if inverse computation was successful
  virtual bool inverse(DataType x_phy, DataType y_phy, DataType z_phy,
                       DataType &x_ref, DataType &y_ref,
                       DataType &z_ref) const = 0;

  /// Given reference coordinates, this routine computes the physical x
  /// coordinates
  virtual DataType x(const Coord &coord_ref) const = 0;
  /// \brief Given reference coordinates, this routine computes the derivatives
  /// in
  ///        direction x of the mapping (ref_coordinates to physical x value)
  virtual DataType x_x(const Coord &coord_ref) const = 0;
  /// \brief Given reference coordinates, this routine computes the derivatives
  /// in
  ///        direction y of the mapping (ref_coordinates to physical x value)
  virtual DataType x_y(const Coord &coord_ref) const = 0;
  /// \brief Given reference coordinates, this routine computes the derivatives
  /// in
  ///        direction z of the mapping (ref_coordinates to physical x value).
  ///        Return value is a dummy value for 2D problems, but neccessary for
  ///        UnitIntegrator

  virtual DataType x_z(const Coord &coord_ref) const { return 0.; }
  /// \brief Given reference coordinates, this routine computes the second
  /// derivatives
  ///        in direction xx of the mapping (ref_coordinates to physical x
  ///        value)
  virtual DataType x_xx(const Coord &coord_ref) const = 0;
  /// \brief Given reference coordinates, this routine computes the second
  /// derivatives
  ///        in direction xy of the mapping (ref_coordinates to physical x
  ///        value)
  virtual DataType x_xy(const Coord &coord_ref) const = 0;
  /// \brief Given reference coordinates, this routine computes the second
  /// derivatives
  ///        in direction xz of the mapping (ref_coordinates to physical x
  ///        value). Return value is a dummy value for 2D problems, but
  ///        neccessary for UnitIntegrator

  virtual DataType x_xz(const Coord &coord_ref) const { return 1.; }
  /// \brief Given reference coordinates, this routine computes the second
  /// derivatives
  ///        in direction yy of the mapping (ref_coordinates to physical x
  ///        value)
  virtual DataType x_yy(const Coord &coord_ref) const = 0;
  /// \brief Given reference coordinates, this routine computes the second
  /// derivatives
  ///        in direction yz of the mapping (ref_coordinates to physical x
  ///        value). Return value is a dummy value for 2D problems, but
  ///        neccessary for UnitIntegrator

  virtual DataType x_yz(const Coord &coord_ref) const { return 1.; }
  /// \brief Given reference coordinates, this routine computes the second
  /// derivatives
  ///        in direction zz of the mapping (ref_coordinates to physical x
  ///        value). Return value is a dummy value for 2D problems, but
  ///        neccessary for UnitIntegrator

  virtual DataType x_zz(const Coord &coord_ref) const { return 1.; }

  /// Given reference coordinates, this computes the physical y coordinates
  virtual DataType y(const Coord &coord_ref) const = 0;
  /// \brief Given reference coordinates, this routine computes the derivatives
  /// in
  ///        direction x of the mapping (ref_coordinates to physical y value)
  virtual DataType y_x(const Coord &coord_ref) const = 0;
  /// \brief Given reference coordinates, this routine computes the derivatives
  /// in
  ///        direction y of the mapping (ref_coordinates to physical y value)
  virtual DataType y_y(const Coord &coord_ref) const = 0;
  /// \brief Given reference coordinates, this routine computes the derivatives
  /// in
  ///        direction z of the mapping (ref_coordinates to physical y value).
  ///        Return value is a dummy value for 2D problems, but neccessary for
  ///        UnitIntegrator

  virtual DataType y_z(const Coord &coord_ref) const { return 0.; }
  /// \brief Given reference coordinates, this routine computes the second
  /// derivatives
  ///        in direction xx of the mapping (ref_coordinates to physical y
  ///        value)
  virtual DataType y_xx(const Coord &coord_ref) const = 0;
  /// \brief Given reference coordinates, this routine computes the second
  /// derivatives
  ///        in direction xy of the mapping (ref_coordinates to physical y
  ///        value)
  virtual DataType y_xy(const Coord &coord_ref) const = 0;
  /// \brief Given reference coordinates, this routine computes the second
  /// derivatives
  ///        in direction xz of the mapping (ref_coordinates to physical y
  ///        value). Return value is a dummy value for 2D problems, but
  ///        neccessary for UnitIntegrator

  virtual DataType y_xz(const Coord &coord_ref) const { return 1.; }
  /// \brief Given reference coordinates, this routine computes the second
  /// derivatives
  ///        in direction yy of the mapping (ref_coordinates to physical y
  ///        value)
  virtual DataType y_yy(const Coord &coord_ref) const = 0;
  /// \brief Given reference coordinates, this routine computes the second
  /// derivatives
  ///        in direction yz of the mapping (ref_coordinates to physical y
  ///        value). Return value is a dummy value for 2D problems, but
  ///        neccessary for UnitIntegrator

  virtual DataType y_yz(const Coord &coord_ref) const { return 1.; }
  /// \brief Given reference coordinates, this routine computes the second
  /// derivatives
  ///        in direction zz of the mapping (ref_coordinates to physical y
  ///        value). Return value is a dummy value for 2D problems, but
  ///        neccessary for UnitIntegrator

  virtual DataType y_zz(const Coord &coord_ref) const { return 1.; }

  /// \brief Given reference coordinates, this computes the physical z
  /// coordinates.
  ///        Return value is a dummy value for 2D problems, but neccessary for
  ///        UnitIntegrator

  virtual DataType z(const Coord &coord_ref) const { return 1.; }
  /// \brief Given reference coordinates, this routine computes the derivatives
  /// in
  ///        direction x of the mapping (ref_coordinates to physical z value).
  ///        Return value is a dummy value for 2D problems, but neccessary for
  ///        UnitIntegrator

  virtual DataType z_x(const Coord &coord_ref) const { return 0.; }
  /// \brief Given reference coordinates, this routine computes the derivatives
  /// in
  ///        direction y of the mapping (ref_coordinates to physical z value).
  ///        Return value is a dummy value for 2D problems, but neccessary for
  ///        UnitIntegrator

  virtual DataType z_y(const Coord &coord_ref) const { return 0.; }
  /// \brief Given reference coordinates, this routine computes the derivatives
  /// in
  ///        direction z of the mapping (ref_coordinates to physical z value).
  ///        Return value is a dummy value for 2D problems, but neccessary for
  ///        UnitIntegrator

  virtual DataType z_z(const Coord &coord_ref) const { return 1.; }
  /// \brief Given reference coordinates, this routine computes the second
  /// derivatives
  ///        in direction xx of the mapping (ref_coordinates to physical z
  ///        value). Return value is a dummy value for 2D problems, but
  ///        neccessary for UnitIntegrator

  virtual DataType z_xx(const Coord &coord_ref) const { return 1.; }
  /// \brief Given reference coordinates, this routine computes the second
  /// derivatives
  ///        in direction xy of the mapping (ref_coordinates to physical z
  ///        value). Return value is a dummy value for 2D problems, but
  ///        neccessary for UnitIntegrator

  virtual DataType z_xy(const Coord &coord_ref) const { return 1.; }
  /// \brief Given reference coordinates, this routine computes the second
  /// derivatives
  ///        in direction xz of the mapping (ref_coordinates to physical z
  ///        value). Return value is a dummy value for 2D problems, but
  ///        neccessary for UnitIntegrator

  virtual DataType z_xz(const Coord &coord_ref) const { return 1.; }
  /// \brief Given reference coordinates, this routine computes the second
  /// derivatives
  ///        in direction yy of the mapping (ref_coordinates to physical z
  ///        value). Return value is a dummy value for 2D problems, but
  ///        neccessary for UnitIntegrator

  virtual DataType z_yy(const Coord &coord_ref) const { return 1.; }
  /// \brief Given reference coordinates, this routine computes the second
  /// derivatives
  ///        in direction yz of the mapping (ref_coordinates to physical z
  ///        value). Return value is a dummy value for 2D problems, but
  ///        neccessary for UnitIntegrator

  virtual DataType z_yz(const Coord &coord_ref) const { return 1.; }
  /// \brief Given reference coordinates, this routine computes the second
  /// derivatives
  ///        in direction zz of the mapping (ref_coordinates to physical z
  ///        value). Return value is a dummy value for 2D problems, but
  ///        neccessary for UnitIntegrator

  virtual DataType z_zz(const Coord &coord_ref) const { return 1.; }

  virtual void transform (const Coord& coord_ref, DataType& x) const {
    x = this->x(coord_ref);
  }
   
  virtual void transform ( const Coord& coord_ref, DataType& x, DataType& y ) const {
    x = this->x ( coord_ref );
    y = this->y ( coord_ref );
  }

  virtual void transform ( const Coord& coord_ref, DataType& x, DataType& y, DataType& z ) const {
    x = this->x ( coord_ref );
    y = this->y ( coord_ref );
    z = this->z ( coord_ref );
  }

  virtual void transform ( const Coord& coord_ref, Coord& coord_mapped ) const {     
    assert ( coord_mapped.size ( ) == coord_ref.size ( ) );
    coord_mapped[0] = this->x ( coord_ref );
    if (this->gdim_ > 1 )
    {
      coord_mapped[1] = this->y ( coord_ref );
      if (this->gdim_ > 2 )
      {
        coord_mapped[2] = this->z ( coord_ref );
      }
    }
   }
            
  /// \brief Check whether a given point is contained in the closure of the
  /// cell.
  ///
  /// \param[in] coord_ref    reference coordinates of the point
  /// \returns  True if reference coordinates are contained in the cell.
  virtual bool contains_reference_point(const Coord &coord_ref) const = 0;

  /// \brief Check whether a given point is contained in the closure of the
  /// cell.
  /// \param[in]  coord_phys   Physical coordinates of the point.
  /// \param[out] coord_ref    Optional output parameter, where reference
  /// coordinates in the cell are stored, if the function returns true. \returns
  /// True if reference coordinates are contained in the cell.
  virtual bool contains_physical_point(const Coord &coord_phys,
                                       Coord *coord_ref) const;

protected:
  /// Vector, which holds the coordinares of every vertex of the physical cell
  Coord coord_vtx_;

  /// Geometrical dimension
  int gdim_;

  /// \details The vector index is calculated by an offset of the
  ///          magnitude i * geometrical dimension
  /// \param[in] i index of vertex id
  /// \param[in] j index of coordinate id (0 for x, 1 for y and 2 for z)
  /// \return index for vector of coordinates coord_vtx_

  inline int ij2ind(int i, int j) const { return i * this->gdim_ + j; }

  /// \brief Computation of inverse mapping in 2d, using Newton's method.
  bool inverse_newton_2d(DataType x_phy, DataType y_phy, DataType &x_ref,
                         DataType &y_ref, DataType x_ref_0 = 0.5,
                         DataType y_ref_0 = 0.5) const;

  /// \brief Computation of inverse mapping in 3d, using Newton's method.
  bool inverse_newton_3d(DataType x_phy, DataType y_phy, DataType z_phy,
                         DataType &x_ref, DataType &y_ref, DataType &z_ref,
                         DataType x_ref_0 = 0.1154, DataType y_ref_0 = 0.1832,
                         DataType z_ref_0 = 0.1385) const;

  /// tolerance for checking whether reference coordinates lie in cell
  DataType eps_;
};

} // namespace doffem
} // namespace hiflow

#endif
