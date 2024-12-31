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

#include "fem/femanager.h"
#include "mesh/iterator.h"

#include "alignedhexahedrontransformation.h"
#include "alignedquadtransformation.h"
#include "bilinearquadtransformation.h"
#include "felagrange.h"
#include "felagrange_hex.h"
#include "felagrange_line.h"
#include "felagrange_pyr.h"
#include "felagrange_quad.h"
#include "felagrange_tet.h"
#include "felagrange_tri.h"
#include "linearlinetransformation.h"
#include "linearpyramidtransformation.h"
#include "lineartetrahedrontransformation.h"
#include "lineartriangletransformation.h"
#include "mesh/periodicity_tools.h"
#include "trilinearhexahedrontransformation.h"

/// \author Michael Schick<br>Martin Baumann<br>Philipp Gerstner

namespace hiflow {
namespace doffem {

template < class DataType >
FEManager< DataType >::FEManager(int dim, int nb_var)
    : dim_(dim), nb_var_(nb_var), num_entities_(-1), mesh_(nullptr) {
  initialized_.resize(nb_var, false);

  // Setting default value true for a continuous ansatz on every variable
  ca_.resize(nb_var_, true);

  fe_tank_.clear();
  fe_inst_.clear();
}

template < class DataType > FEManager< DataType >::~FEManager() {
  for (size_t i = 0, e_i = cell_transformation_.size(); i != e_i; ++i) {
    delete cell_transformation_[i];
  }
}

template < class DataType >
void FEManager< DataType >::set_mesh(const mesh::Mesh &mesh) {
  mesh_ = &mesh;
  num_entities_ = mesh_->num_entities(dim_);
  init_cell_transformation();
  fe_tank_.resize(nb_var_ * num_entities_);
}

/// \details By iterating through the given mesh, for each mesh cell the
/// coordinates
///          are extracted and an instance of a transformation with correct type
///          is stored.

template < class DataType >
void FEManager< DataType >::init_cell_transformation() {
  assert(mesh_ != 0);

  cell_transformation_.clear();
  cell_transformation_.resize(num_entities_);

  std::vector< mesh::MasterSlave > period = mesh_->get_period();

  for (mesh::EntityIterator it = mesh_->begin(dim_), e_it = mesh_->end(dim_);
       it != e_it; ++it) {
    // Get Cell Type (line, hex, quad, tri, tet....)

    mesh::CellType::Tag cell_type = it->cell_type().tag();
    mesh::AlignNumber align_number = it->get_align_number();

    // Initialize Cell Transformations

    switch (cell_type) {

    case mesh::CellType::LINE: {
      cell_transformation_[it->index()] =
          new LinearLineTransformation< DataType >(1);
      Coord coord_vtx;
      it->get_coordinates(coord_vtx);

      if (!period.empty()) {
        Coord tmp_coord = coord_vtx;
        coord_vtx = unperiodify(tmp_coord, 1, period);
      }

      cell_transformation_[it->index()]->reinit(coord_vtx);
      break;
    }

    case mesh::CellType::TRIANGLE: {
      cell_transformation_[it->index()] =
          new LinearTriangleTransformation< DataType >(2);
      Coord coord_vtx;
      it->get_coordinates(coord_vtx);

      if (!period.empty()) {
        Coord tmp_coord = coord_vtx;
        coord_vtx = unperiodify(tmp_coord, 2, period);
      }

      cell_transformation_[it->index()]->reinit(coord_vtx);
      break;
    }

    case mesh::CellType::QUADRILATERAL: {
      // TODO: differentiate more align cases
      if (align_number >= 4) {
        // 4: aligned w.r.t x and y axis
        cell_transformation_[it->index()] =
            new AlignedQuadTransformation< DataType >(2);
      } else {
        cell_transformation_[it->index()] =
            new BiLinearQuadTransformation< DataType >(2);
      }
      Coord coord_vtx;
      it->get_coordinates(coord_vtx);

      if (!period.empty()) {
        Coord tmp_coord = coord_vtx;
        coord_vtx = unperiodify(tmp_coord, 2, period);
      }

      cell_transformation_[it->index()]->reinit(coord_vtx);
      break;
    }

    case mesh::CellType::TETRAHEDRON: {
      cell_transformation_[it->index()] =
          new LinearTetrahedronTransformation< DataType >(3);
      Coord coord_vtx;
      it->get_coordinates(coord_vtx);

      if (!period.empty()) {
        Coord tmp_coord = coord_vtx;
        coord_vtx = unperiodify(tmp_coord, 3, period);
      }

      cell_transformation_[it->index()]->reinit(coord_vtx);
      break;
    }

    case mesh::CellType::HEXAHEDRON: {
      // TODO: differentiate more align cases
      if (align_number >= 7) {
        // 7: aligned w.r.t x, y and z axis
        cell_transformation_[it->index()] =
            new AlignedHexahedronTransformation< DataType >(3);
      } else {
        cell_transformation_[it->index()] =
            new TriLinearHexahedronTransformation< DataType >(3);
      }
      Coord coord_vtx;
      it->get_coordinates(coord_vtx);

      if (!period.empty()) {
        Coord tmp_coord = coord_vtx;
        coord_vtx = unperiodify(tmp_coord, 3, period);
      }

      cell_transformation_[it->index()]->reinit(coord_vtx);
      break;
    }

    case mesh::CellType::PYRAMID: {
      cell_transformation_[it->index()] =
          new LinearPyramidTransformation< DataType >(3);
      Coord coord_vtx;
      it->get_coordinates(coord_vtx);

      if (!period.empty()) {
        Coord tmp_coord = coord_vtx;
        coord_vtx = unperiodify(tmp_coord, 3, period);
      }

      cell_transformation_[it->index()]->reinit(coord_vtx);
      break;
    }

    default:
      assert(0);
    }
  }
}

/// \details Initializing is done by a loop over all variables

template < class DataType >
void FEManager< DataType >::init_fe_tank(
    const typename FEType< DataType >::FEAnsatz &ansatz,
    const std::vector< std::vector< int > > &param) {
  assert(param.size() == nb_var_);
  for (int var = 0; var < nb_var_; ++var) {
    init_fe_tank(var, ansatz, param[var]);
  }
}

/// \details By iterating through the given mesh, on each cell for a given
/// variable the
///          desired finite elements are referenced by FEInstance \see
///          FEInstance. The protected variable fe_tank_ is filled with pointers
///          to the corresponding ansatz.

template < class DataType >
void FEManager< DataType >::init_fe_tank(
    int var, const typename FEType< DataType >::FEAnsatz &ansatz,
    const std::vector< int > &param) {
  assert(!initialized_[var]);

  assert(!param.empty());

  // Initialize pointers to existing entity types
  mesh::EntityNumber line_representer = -1;
  mesh::EntityNumber tri_representer = -1;
  mesh::EntityNumber quad_representer = -1;
  mesh::EntityNumber tet_representer = -1;
  mesh::EntityNumber hex_representer = -1;
  mesh::EntityNumber pyr_representer = -1;

  const size_t offset = var * num_entities_;

  for (mesh::EntityIterator it = mesh_->begin(dim_), e_it = mesh_->end(dim_);
       it != e_it; ++it) {
    // 1st: Get Cell Type (Tet, Hex, Quad ...)
    mesh::CellType::Tag cell_type = it->cell_type().tag();

    // 2nd: Store Finite Element in FE Tank
    switch (cell_type) {

    case mesh::CellType::LINE: {
      if (ansatz == FEType< DataType >::LAGRANGE) {
        line_representer = it->index();
        fe_tank_[it->index() + offset] =
            fe_inst_.get_lagrange_line_element(param[0]);
      } else {
        assert(0);
      }
      break;
    }

    case mesh::CellType::TRIANGLE: {
      if (ansatz == FEType< DataType >::LAGRANGE) {
        tri_representer = it->index();
        fe_tank_[it->index() + offset] =
            fe_inst_.get_lagrange_tri_element(param[0]);
      } else {
        assert(0);
      }
      break;
    }

    case mesh::CellType::QUADRILATERAL: {
      if (ansatz == FEType< DataType >::LAGRANGE) {
        quad_representer = it->index();
        fe_tank_[it->index() + offset] =
            fe_inst_.get_lagrange_quad_element(param[0]);
      } else {
        assert(0);
      }
      break;
    }

    case mesh::CellType::TETRAHEDRON: {
      if (ansatz == FEType< DataType >::LAGRANGE) {
        tet_representer = it->index();
        fe_tank_[it->index() + offset] =
            fe_inst_.get_lagrange_tet_element(param[0]);
      } else
        assert(0);
      break;
    }

    case mesh::CellType::HEXAHEDRON: {
      if (ansatz == FEType< DataType >::LAGRANGE) {
        hex_representer = it->index();
        fe_tank_[it->index() + offset] =
            fe_inst_.get_lagrange_hex_element(param[0]);
      } else {
        assert(0);
      }
      break;
    }

    case mesh::CellType::PYRAMID: {
      if (ansatz == FEType< DataType >::LAGRANGE) {
        pyr_representer = it->index();
        fe_tank_[it->index() + offset] =
            fe_inst_.get_lagrange_pyr_element(param[0]);
      } else {
        assert(0);
      }
      break;
    }

    default:
      assert(0);
    }
  }

  // At least one cell should have been found
  assert(!(line_representer == -1 && tri_representer == -1 &&
           quad_representer == -1 && tet_representer == -1 &&
           hex_representer == -1 && pyr_representer == -1));

  // Last:  setting up numbering of reference cell via transformation of an
  // arbitrary
  //        (here the last) mesh entity back to reference cell and iteration
  //        through the vertex ids

  if (line_representer != -1) {
    std::vector< DataType > coordinates;
    mesh::Entity entity = mesh_->get_entity(dim_, line_representer);
    entity.get_coordinates(coordinates);

    LinearLineTransformation< DataType > llt(1);
    // std::vector<Coord> ref_coords = llt.get_reference_coordinates();

    llt.reinit(coordinates);

    std::vector< Coord > ref_coords(2);
    for (size_t i = 0; i < 2; ++i) {
      ref_coords[i].resize(1);
    }

    for (size_t i = 0; i < 2; ++i) {
      bool found = llt.inverse(coordinates[i], ref_coords[i][0]);
      assert(found);
    }

    fe_inst_.init_line_elements(ref_coords);
  }

  if (tri_representer != -1) {
    std::vector< DataType > coordinates;
    mesh::Entity entity = mesh_->get_entity(dim_, tri_representer);
    entity.get_coordinates(coordinates);

    LinearTriangleTransformation< DataType > ltt(2);
    // std::vector<Coord> ref_coords = ltt.get_reference_coordinates();

    ltt.reinit(coordinates);

    std::vector< Coord > ref_coords(3);
    for (size_t i = 0; i < 3; ++i) {
      ref_coords[i].resize(2);
    }

    for (size_t i = 0; i < 3; ++i) {
      bool found = ltt.inverse(coordinates[2 * i], coordinates[2 * i + 1],
                               ref_coords[i][0], ref_coords[i][1]);
      assert(found);
    }

    fe_inst_.init_triangle_elements(ref_coords);
  }

  if (quad_representer != -1) {
    std::vector< DataType > coordinates;
    mesh::Entity entity = mesh_->get_entity(dim_, quad_representer);
    entity.get_coordinates(coordinates);

    std::vector< Coord > ref_coords(4);
    for (size_t i = 0; i < 4; ++i) {
      ref_coords[i].resize(2);
    }

    // TODO: what if not all cells are axis aligned, i.e if we have both
    // BiLinear and Aligned Quad trafos?
    if (!this->mesh_->is_rectangular()) {
      BiLinearQuadTransformation< DataType > blqt(2);
      blqt.reinit(coordinates);
      for (size_t i = 0; i < 4; ++i) {
        bool found = blqt.inverse(coordinates[2 * i], coordinates[2 * i + 1],
                                  ref_coords[i][0], ref_coords[i][1]);
        assert(found);
      }
    } else {
      AlignedQuadTransformation< DataType > blqt(2);
      blqt.reinit(coordinates);
      for (size_t i = 0; i < 4; ++i) {
        bool found = blqt.inverse(coordinates[2 * i], coordinates[2 * i + 1],
                                  ref_coords[i][0], ref_coords[i][1]);
        assert(found);
      }
    }

    fe_inst_.init_quadrilateral_elements(ref_coords);
  }

  if (tet_representer != -1) {
    std::vector< DataType > coordinates;
    mesh::Entity entity = mesh_->get_entity(dim_, tet_representer);
    entity.get_coordinates(coordinates);

    LinearTetrahedronTransformation< DataType > ltt(3);

    ltt.reinit(coordinates);

    std::vector< Coord > ref_coords(4);
    for (size_t i = 0; i < 4; ++i) {
      ref_coords[i].resize(3);
    }

    for (size_t i = 0; i < 4; ++i) {
      bool found = ltt.inverse(coordinates[3 * i], coordinates[3 * i + 1],
                               coordinates[3 * i + 2], ref_coords[i][0],
                               ref_coords[i][1], ref_coords[i][2]);
      assert(found);
    }

    fe_inst_.init_tetrahedron_elements(ref_coords);
  }

  if (hex_representer != -1) {
    std::vector< DataType > coordinates;
    mesh::Entity entity = mesh_->get_entity(dim_, hex_representer);
    entity.get_coordinates(coordinates);

    std::vector< Coord > ref_coords(8);
    for (size_t i = 0; i < 8; ++i) {
      ref_coords[i].resize(3);
    }

    // TODO: what if not all cells are axis aligned, i.e if we have both
    // TriLinear and Aligned Hexahedon trafos?
    if (!this->mesh_->is_rectangular()) {
      TriLinearHexahedronTransformation< DataType > tlht(3);
      tlht.reinit(coordinates);
      for (size_t i = 0; i < 8; ++i) {
        bool found = tlht.inverse(coordinates[3 * i], coordinates[3 * i + 1],
                                  coordinates[3 * i + 2], ref_coords[i][0],
                                  ref_coords[i][1], ref_coords[i][2]);
        assert(found);
      }
    } else {
      AlignedHexahedronTransformation< DataType > tlht(3);
      tlht.reinit(coordinates);
      for (size_t i = 0; i < 8; ++i) {
        bool found = tlht.inverse(coordinates[3 * i], coordinates[3 * i + 1],
                                  coordinates[3 * i + 2], ref_coords[i][0],
                                  ref_coords[i][1], ref_coords[i][2]);
        assert(found);
      }
    }

    fe_inst_.init_hexahedron_elements(ref_coords);
  }

  if (pyr_representer != -1) {
    std::vector< DataType > coordinates;
    mesh::Entity entity = mesh_->get_entity(dim_, pyr_representer);
    entity.get_coordinates(coordinates);

    LinearPyramidTransformation< DataType > lpr(3);
    // std::vector<Coord> ref_coords = lpr.get_reference_coordinates();

    lpr.reinit(coordinates);

    std::vector< Coord > ref_coords(5);
    for (size_t i = 0; i < 5; ++i) {
      ref_coords[i].resize(3);
    }

    for (size_t i = 0; i < 5; ++i) {
      bool found = lpr.inverse(coordinates[3 * i], coordinates[3 * i + 1],
                               coordinates[3 * i + 2], ref_coords[i][0],
                               ref_coords[i][1], ref_coords[i][2]);
      assert(found);
    }

    fe_inst_.init_pyramid_elements(ref_coords);
  }

  initialized_[var] = true;
}

template < class DataType > void FEManager< DataType >::get_status() const {
  std::cout << "DIM:   " << dim_ << std::endl;
  std::cout << "#vars: " << nb_var_ << std::endl;

  fe_inst_.get_status();

  std::cout << "tank size: " << fe_tank_.size() << std::endl;
  for (size_t i = 0; i < fe_tank_.size(); ++i) {
    std::cout << "\t" << i << "\t" << fe_tank_[i]->get_name() << std::endl;
  }
}

template class FEManager< double >;
template class FEManager< float >;
} // namespace doffem
} // namespace hiflow
