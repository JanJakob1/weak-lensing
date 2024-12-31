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

#ifndef __FEM_FEINSTANCE_H_
#define __FEM_FEINSTANCE_H_

#include "dof/dof_fem_types.h"
#include <vector>
#include "common/vector_algebra.h"
#include "fem/felagrange_hex.h"
#include "fem/felagrange_line.h"
#include "fem/felagrange_pyr.h"
#include "fem/felagrange_quad.h"
#include "fem/felagrange_tet.h"
#include "fem/felagrange_tri.h"

namespace hiflow {
namespace doffem {

template < class DataType, int DIM > class FELagrangeLine;

template < class DataType, int DIM > class FELagrangeTri;

template < class DataType, int DIM > class FELagrangeQuad;

template < class DataType, int DIM > class FELagrangeTet;

template < class DataType, int DIM > class FELagrangeHex;

template < class DataType, int DIM > class FELagrangePyr;

///
/// \class FEInstance feinstance.h
/// \brief Holds instances of different FE ansatz functions
/// \author Michael Schick<br>Martin Baumann
///

template < class DataType, int DIM > class FEInstance {
public:
  typedef Vec<DIM, DataType> Coord;

  /// Default Constructor
  FEInstance();
  /// Default Destructor
  ~FEInstance();

  /// Returns pointer to an instance of a Lagrangian Line element of given
  /// degree
  FELagrangeLine< DataType, DIM > *get_lagrange_line_element(int degree);
  /// Returns pointer to an instance of a Lagrangian Triangle element of given
  /// degree
  FELagrangeTri< DataType, DIM > *get_lagrange_tri_element(int degree);
  /// Returns pointer to an instance of a Lagrangian Quadrilateral element of
  /// given degree
  FELagrangeQuad< DataType, DIM > *get_lagrange_quad_element(int degree);
  /// Returns pointer to an instance of a Lagrangian Tetrahedron element of
  /// given degree
  FELagrangeTet< DataType, DIM > *get_lagrange_tet_element(int degree);
  /// Returns pointer to an instance of a Lagrangian Hexahedron element of given
  /// degree
  FELagrangeHex< DataType, DIM > *get_lagrange_hex_element(int degree);
  /// Returns pointer to an instance of a Lagrangian Pyramid element of given
  /// degree
  FELagrangePyr< DataType, DIM > *get_lagrange_pyr_element(int degree);

  /// Returns pointer to an instance of a Lagrangian Line element of given
  /// degree (const)
  FELagrangeLine< DataType, DIM > *get_const_lagrange_line_element(int degree) const;
  /// Returns pointer to an instance of a Lagrangian Triangle element of given
  /// degree (const)
  FELagrangeTri< DataType, DIM > *get_const_lagrange_tri_element(int degree) const;
  /// Returns pointer to an instance of a Lagrangian Quadrilateral element of
  /// given degree (const)
  FELagrangeQuad< DataType, DIM > *get_const_lagrange_quad_element(int degree) const;
  /// Returns pointer to an instance of a Lagrangian Tetrahedron element of
  /// given degree (const)
  FELagrangeTet< DataType, DIM > *get_const_lagrange_tet_element(int degree) const;
  /// Returns pointer to an instance of a Lagrangian Hexahedron element of given
  /// degree (const)
  FELagrangeHex< DataType, DIM > *get_const_lagrange_hex_element(int degree) const;
  /// Returns pointer to an instance of a Lagrangian Pyramid element of given
  /// degree (const)
  FELagrangePyr< DataType, DIM > *get_const_lagrange_pyr_element(int degree) const;

  /// Return status information about FEInstance
  void get_status() const;
  /// Return number of initialized Finite Elements (use only after complete
  /// initialization!)

  int nfe() const { return id_; }

  /// Clearing all stored finite elements and reseting FEInstance
  void clear();

  /// Initialize the coordinates on the reference line for vertices
  void init_line_elements(const std::vector< Coord > &ref_coords);
  /// Initialize the coordinates on the reference triangle for vertices
  void init_triangle_elements(const std::vector< Coord > &ref_coords);
  /// Initialize the coordinates on the reference quadrilateral for vertices
  void init_quadrilateral_elements(const std::vector< Coord > &ref_coords);
  /// Initialize the coordinates on the reference tetrahedron for vertices
  void init_tetrahedron_elements(const std::vector< Coord > &ref_coords);
  /// Initialize the coordinates on the reference hexahedron for vertices
  void init_hexahedron_elements(const std::vector< Coord > &ref_coords);
  /// Initialize the coordinates on the reference pyramid for vertices
  void init_pyramid_elements(const std::vector< Coord > &ref_coords);

  void pass_coord_to_element ( FELagrangeLine<DataType, DIM >* fe_lag, 
                               const std::vector< Coord > &ref_coords );
  void pass_coord_to_element ( FELagrangeTri<DataType, DIM >* fe_lag, 
                               const std::vector< Coord > &ref_coords );
  void pass_coord_to_element ( FELagrangeQuad<DataType, DIM >* fe_lag, 
                               const std::vector< Coord > &ref_coords );
  void pass_coord_to_element ( FELagrangeTet<DataType, DIM >* fe_lag, 
                               const std::vector< Coord > &ref_coords );
  void pass_coord_to_element ( FELagrangeHex<DataType, DIM >* fe_lag, 
                               const std::vector< Coord > &ref_coords );
  void pass_coord_to_element ( FELagrangePyr<DataType, DIM >* fe_lag, 
                               const std::vector< Coord > &ref_coords );

private:
  /// Id counter
  int id_;

  /// Stored Lagrange instances on a Line
  std::vector< FELagrangeLine< DataType, DIM > * > fe_lag_line_;

  /// Stored Lagrange instances on a Triangle
  std::vector< FELagrangeTri< DataType, DIM > * > fe_lag_tri_;

  /// Stored Lagrange instances on a Quadrilateral
  std::vector< FELagrangeQuad< DataType, DIM > * > fe_lag_quad_;

  /// Stored Lagrange instances on a Tetrahedron
  std::vector< FELagrangeTet< DataType, DIM > * > fe_lag_tet_;

  /// Stored Lagrange instances on a Hexahedron
  std::vector< FELagrangeHex< DataType, DIM > * > fe_lag_hex_;

  /// Stored Lagrange instances on a Pyramid
  std::vector< FELagrangePyr< DataType, DIM > * > fe_lag_pyr_;
};

template < class DataType, int DIM > FEInstance< DataType, DIM >::FEInstance() { id_ = 0; }

template < class DataType, int DIM > FEInstance< DataType, DIM >::~FEInstance() { clear(); }

/// \details Every allocated instance of a finite element is being
///          deleted and the vectors storing this information are cleared

template < class DataType, int DIM > void FEInstance< DataType, DIM >::clear() {
  for (size_t i = 0, e_i = fe_lag_line_.size(); i != e_i; ++i) {
    delete fe_lag_line_[i];
  }

  for (size_t i = 0, e_i = fe_lag_tri_.size(); i != e_i; ++i) {
    delete fe_lag_tri_[i];
  }

  for (size_t i = 0, e_i = fe_lag_quad_.size(); i != e_i; ++i) {
    delete fe_lag_quad_[i];
  }

  for (size_t i = 0, e_i = fe_lag_tet_.size(); i != e_i; ++i) {
    delete fe_lag_tet_[i];
  }

  for (size_t i = 0, e_i = fe_lag_hex_.size(); i != e_i; ++i) {
    delete fe_lag_hex_[i];
  }

  for (size_t i = 0, e_i = fe_lag_pyr_.size(); i != e_i; ++i) {
    delete fe_lag_pyr_[i];
  }

  fe_lag_line_.clear();
  fe_lag_tri_.clear();
  fe_lag_quad_.clear();
  fe_lag_tet_.clear();
  fe_lag_hex_.clear();
  fe_lag_pyr_.clear();

  id_ = 0;
}

template < class DataType, int DIM >
void FEInstance< DataType, DIM >::init_line_elements(const std::vector< Coord > &ref_coords) {
  assert(ref_coords.size() == 2);

  for (size_t deg = 0, e_deg = fe_lag_line_.size(); deg != e_deg; ++deg) {
    if (!fe_lag_line_[deg]->get_init_status()) {
      this->pass_coord_to_element (fe_lag_line_[deg], ref_coords);
      fe_lag_line_[deg]->init();
      fe_lag_line_[deg]->set_id(id_);
      ++id_;
    }
  }
}

template < class DataType, int DIM >
void FEInstance< DataType, DIM >::init_triangle_elements(const std::vector< Coord > &ref_coords) {
  assert(ref_coords.size() == 3);

  for (size_t deg = 0, e_deg = fe_lag_tri_.size(); deg != e_deg; ++deg) {
    if (!fe_lag_tri_[deg]->get_init_status()) {
      this->pass_coord_to_element (fe_lag_tri_[deg], ref_coords);
      fe_lag_tri_[deg]->init();
      fe_lag_tri_[deg]->set_id(id_);
      ++id_;
    }
  }
}

template < class DataType, int DIM >
void FEInstance< DataType, DIM >::init_quadrilateral_elements(const std::vector< Coord > &ref_coords) {
  assert(ref_coords.size() == 4);

  for (size_t deg = 0, e_deg = fe_lag_quad_.size(); deg != e_deg; ++deg) {
    if (!fe_lag_quad_[deg]->get_init_status()) {
      this->pass_coord_to_element (fe_lag_quad_[deg], ref_coords);
      fe_lag_quad_[deg]->init();
      fe_lag_quad_[deg]->set_id(id_);
      ++id_;
    }
  }
}

template < class DataType, int DIM >
void FEInstance< DataType, DIM >::init_tetrahedron_elements(const std::vector< Coord > &ref_coords) {
  assert(ref_coords.size() == 4);

  for (size_t deg = 0, e_deg = fe_lag_tet_.size(); deg != e_deg; ++deg) {
    if (!fe_lag_tet_[deg]->get_init_status()) {
      this->pass_coord_to_element (fe_lag_tet_[deg], ref_coords);
      fe_lag_tet_[deg]->init();
      fe_lag_tet_[deg]->set_id(id_);
      ++id_;
    }
  }
}

template < class DataType, int DIM >
void FEInstance< DataType, DIM >::init_hexahedron_elements(const std::vector< Coord > &ref_coords) {
  assert(ref_coords.size() == 8);

  for (size_t deg = 0, e_deg = fe_lag_hex_.size(); deg != e_deg; ++deg) {
    if (!fe_lag_hex_[deg]->get_init_status()) {
      this->pass_coord_to_element (fe_lag_hex_[deg], ref_coords);
      fe_lag_hex_[deg]->init();
      fe_lag_hex_[deg]->set_id(id_);
      ++id_;
    }
  }
}

template < class DataType, int DIM >
void FEInstance< DataType, DIM >::init_pyramid_elements(const std::vector< Coord > &ref_coords) {
  assert(ref_coords.size() == 5);

  for (size_t deg = 0, e_deg = fe_lag_pyr_.size(); deg != e_deg; ++deg) {
    if (!fe_lag_pyr_[deg]->get_init_status()) {
      this->pass_coord_to_element (fe_lag_pyr_[deg], ref_coords);
      fe_lag_pyr_[deg]->init();
      fe_lag_pyr_[deg]->set_id(id_);
      ++id_;
    }
  }
}

/// \details Checks wether a lagrangian finite element on a line with
///          given degree is already existing. If not, then all elements up to
///          the given degree are created.

template < class DataType, int DIM >
FELagrangeLine< DataType, DIM > *
FEInstance< DataType, DIM >::get_lagrange_line_element(int degree) {
  if (fe_lag_line_.size() <= degree) {
    // add all FETypes up to degree
    for (int i = fe_lag_line_.size(); i < degree + 1; ++i) {
      FELagrangeLine< DataType, DIM > *dummy_lag = new FELagrangeLine< DataType, DIM >;
      dummy_lag->set_my_deg(i);
      fe_lag_line_.push_back(dummy_lag);
    }
  }

  return fe_lag_line_[degree];
}

/// \details Checks wether a lagrangian finite element on a triangle with
///          given degree is already existing. If not, then all elements up to
///          the given degree are created.

template < class DataType, int DIM >
FELagrangeTri< DataType, DIM > *
FEInstance< DataType, DIM >::get_lagrange_tri_element(int degree) {
  if (fe_lag_tri_.size() <= degree) {
    // add all FETypes up to degree
    for (int i = fe_lag_tri_.size(); i < degree + 1; ++i) {
      FELagrangeTri< DataType, DIM > *dummy_lag = new FELagrangeTri< DataType, DIM >;
      dummy_lag->set_my_deg(i);
      fe_lag_tri_.push_back(dummy_lag);
    }
  }

  return fe_lag_tri_[degree];
}

/// \details Checks wether a lagrangian finite element on a quadrilateral with
///          given degree is already existing. If not, then all elements up to
///          the given degree are created.

template < class DataType, int DIM >
FELagrangeQuad< DataType, DIM > *
FEInstance< DataType, DIM >::get_lagrange_quad_element(int degree) {
  if (fe_lag_quad_.size() <= degree) {
    // add all FETypes up to degree
    for (int i = fe_lag_quad_.size(); i < degree + 1; ++i) {
      FELagrangeQuad< DataType, DIM > *dummy_lag = new FELagrangeQuad< DataType, DIM >;
      dummy_lag->set_my_deg(i);
      fe_lag_quad_.push_back(dummy_lag);
    }
  }

  return fe_lag_quad_[degree];
}

/// \details Checks wether a lagrangian finite element on a tetrahedron with
///          given degree is already existing. If not, then all elements up to
///          the given degree are created.

template < class DataType, int DIM >
FELagrangeTet< DataType, DIM > *
FEInstance< DataType, DIM >::get_lagrange_tet_element(int degree) {
  if (fe_lag_tet_.size() <= degree) {
    // add all FETypes up to degree
    for (int i = fe_lag_tet_.size(); i < degree + 1; ++i) {
      FELagrangeTet< DataType, DIM > *dummy_lag = new FELagrangeTet< DataType, DIM >;
      dummy_lag->set_my_deg(i);
      fe_lag_tet_.push_back(dummy_lag);
    }
  }

  return fe_lag_tet_[degree];
}

/// \details Checks wether a lagrangian finite element on a hexahedron with
///          given degree is already existing. If not, then all elements up to
///          the given degree are created.

template < class DataType, int DIM >
FELagrangeHex< DataType, DIM > *
FEInstance< DataType, DIM >::get_lagrange_hex_element(int degree) {
  if (fe_lag_hex_.size() <= degree) {
    // add all FETypes up to degree
    for (int i = fe_lag_hex_.size(); i < degree + 1; ++i) {
      FELagrangeHex< DataType, DIM > *dummy_lag = new FELagrangeHex< DataType, DIM >;
      dummy_lag->set_my_deg(i);
      fe_lag_hex_.push_back(dummy_lag);
    }
  }

  return fe_lag_hex_[degree];
}

/// \details Checks wether a lagrangian finite element on a pyramid with
///          given degree is already existing. If not, then all elements up to
///          the given degree are created.

template < class DataType, int DIM >
FELagrangePyr< DataType, DIM > *
FEInstance< DataType, DIM >::get_lagrange_pyr_element(int degree) {
  if (fe_lag_pyr_.size() <= degree) {
    // add all FETypes up to degree
    for (int i = fe_lag_pyr_.size(); i < degree + 1; ++i) {
      FELagrangePyr< DataType, DIM > *dummy_lag = new FELagrangePyr< DataType, DIM >;
      dummy_lag->set_my_deg(i);
      fe_lag_pyr_.push_back(dummy_lag);
    }
  }

  return fe_lag_pyr_[degree];
}

template < class DataType, int DIM >
FELagrangeLine< DataType, DIM > *
FEInstance< DataType, DIM >::get_const_lagrange_line_element(int degree) const {
  return fe_lag_line_[degree];
}

template < class DataType, int DIM >
FELagrangeTri< DataType, DIM > *
FEInstance< DataType, DIM >::get_const_lagrange_tri_element(int degree) const {
  return fe_lag_tri_[degree];
}

template < class DataType, int DIM >
FELagrangeQuad< DataType, DIM > *
FEInstance< DataType, DIM >::get_const_lagrange_quad_element(int degree) const {
  return fe_lag_quad_[degree];
}

template < class DataType, int DIM >
FELagrangeTet< DataType, DIM > *
FEInstance< DataType, DIM >::get_const_lagrange_tet_element(int degree) const {
  return fe_lag_tet_[degree];
}

template < class DataType, int DIM >
FELagrangeHex< DataType, DIM > *
FEInstance< DataType, DIM >::get_const_lagrange_hex_element(int degree) const {
  return fe_lag_hex_[degree];
}

template < class DataType, int DIM >
FELagrangePyr< DataType, DIM > *
FEInstance< DataType, DIM >::get_const_lagrange_pyr_element(int degree) const {
  return fe_lag_pyr_[degree];
}

template < class DataType, int DIM > 
void FEInstance< DataType, DIM >::get_status() const {
  // std::cout << "FEInstance size: " << fe_lag_hex_.size() << std::endl;

  std::cout << "  Lines" << std::endl;
  for (size_t i = 0; i < fe_lag_line_.size(); ++i) {
    std::cout << "\t" << i << "\t" << fe_lag_line_[i]->get_name() << " ("
              << fe_lag_line_[i]->get_my_deg() << ")"
              << "\t@ " << get_const_lagrange_line_element(i) << std::endl;
  }

  std::cout << "  Triangles" << std::endl;
  for (size_t i = 0; i < fe_lag_tri_.size(); ++i) {
    std::cout << "\t" << i << "\t" << fe_lag_tri_[i]->get_name() << " ("
              << fe_lag_tri_[i]->get_my_deg() << ")"
              << "\t@ " << get_const_lagrange_tri_element(i) << std::endl;
  }

  std::cout << "  Quadrilaterals" << std::endl;
  for (size_t i = 0; i < fe_lag_quad_.size(); ++i) {
    std::cout << "\t" << i << "\t" << fe_lag_quad_[i]->get_name() << " ("
              << fe_lag_quad_[i]->get_my_deg() << ")"
              << "\t@ " << get_const_lagrange_quad_element(i) << std::endl;
  }

  std::cout << "  Tetrahedrons" << std::endl;
  for (size_t i = 0; i < fe_lag_tet_.size(); ++i) {
    std::cout << "\t" << i << "\t" << fe_lag_tet_[i]->get_name() << " ("
              << fe_lag_tet_[i]->get_my_deg() << ")"
              << "\t@ " << get_const_lagrange_tet_element(i) << std::endl;
  }

  std::cout << "  Hexahedrons" << std::endl;
  for (size_t i = 0; i < fe_lag_hex_.size(); ++i) {
    std::cout << "\t" << i << "\t" << fe_lag_hex_[i]->get_name() << " ("
              << fe_lag_hex_[i]->get_my_deg() << ")"
              << "\t@ " << get_const_lagrange_hex_element(i) << std::endl;
  }

  std::cout << "  Pyramid" << std::endl;
  for (size_t i = 0; i < fe_lag_pyr_.size(); ++i) {
    std::cout << "\t" << i << "\t" << fe_lag_pyr_[i]->get_name() << " ("
              << fe_lag_pyr_[i]->get_my_deg() << ")"
              << "\t@ " << get_const_lagrange_pyr_element(i) << std::endl;
  }
}

// line element
template < class DataType, int DIM > 
void FEInstance< DataType, DIM >::pass_coord_to_element ( FELagrangeLine< DataType, DIM >* fe_lag, const std::vector< Vec<DIM,DataType> >&ref_coords ) 
{
  fe_lag->set_ref_vtx_coords(ref_coords);
}

// triangle element
template < class DataType, int DIM > 
void FEInstance< DataType, DIM >::pass_coord_to_element ( FELagrangeTri< DataType, DIM >* fe_lag, const std::vector< Vec<DIM, DataType> >&ref_coords ) 
{
  fe_lag->set_ref_vtx_coords(ref_coords);
}

// quad element
template < class DataType, int DIM > 
void FEInstance< DataType, DIM >::pass_coord_to_element ( FELagrangeQuad< DataType, DIM >* fe_lag, const std::vector< Vec<DIM, DataType> >&ref_coords ) 
{
  fe_lag->set_ref_vtx_coords(ref_coords);
}

// tet element
template < class DataType, int DIM > 
void FEInstance< DataType, DIM >::pass_coord_to_element ( FELagrangeTet< DataType, DIM >* fe_lag, const std::vector< Vec<DIM, DataType> >&ref_coords ) 
{
  fe_lag->set_ref_vtx_coords(ref_coords);
}

// hex element
template < class DataType, int DIM > 
void FEInstance< DataType, DIM >::pass_coord_to_element ( FELagrangeHex< DataType, DIM >* fe_lag, const std::vector< Vec<DIM, DataType> >&ref_coords ) 
{
  fe_lag->set_ref_vtx_coords(ref_coords);
}

// pyr element
template < class DataType, int DIM > 
void FEInstance< DataType, DIM >::pass_coord_to_element ( FELagrangePyr< DataType, DIM >* fe_lag, const std::vector< Vec<DIM, DataType> >&ref_coords ) 
{
  fe_lag->set_ref_vtx_coords(ref_coords);
}

} // namespace doffem
} // namespace hiflow
#endif
