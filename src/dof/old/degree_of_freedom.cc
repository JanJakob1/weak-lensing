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

#include "degree_of_freedom.h"

#include <mpi.h>

#include "common/log.h"
#include "common/macros.h"
#include "fem/femanager.h"
#include "mesh/iterator.h"

#include "numbering_lagrange.h"

/// \author Michael Schick<br>Martin Baumann

namespace hiflow {
namespace doffem {

template < class DataType > bool DegreeOfFreedom< DataType >::check_mesh() {
  bool ret = true;

  // loop over mesh cells
  for (mesh::EntityIterator it = this->mesh_->begin(this->get_dim()),
       e_it = this->mesh_->end(this->get_dim());
       it != e_it; ++it) {
    // get coordinates of current cell
    std::vector< DataType > coords;
    it->get_coordinates(coords);
    // check whether vertex coordinates are conform to numbering rule of cell
    // type
    if (!static_cast< bool >(
          it->cell_type().check_cell_geometry(coords, it->gdim()))) {
      std::cout << "ERROR: vertex numbering of cell " << it->index()
                << " is not conform to cell type." << std::endl;
      ret = false;
    }
  }
  return ret;
}

template < class DataType > DegreeOfFreedom< DataType >::DegreeOfFreedom() {
  this->mesh_ = nullptr;

  this->tdim_ = -1;
  this->nvar_ = -1;

  this->number_of_dofs_total_ = -1;

  this->mesh_flag_ = false;
  this->fe_manager_flag_ = false;
  this->strategy_ = false;

  this->number_strategy_ = nullptr;
}

template < class DataType > DegreeOfFreedom< DataType >::~DegreeOfFreedom() {
  delete this->number_strategy_;
}

template < class DataType >
void DegreeOfFreedom< DataType >::init_numbering_strategy(
  const std::string &number_strategy) {
  assert(this->mesh_flag_);
  assert(this->fe_manager_flag_);

  if (number_strategy.empty() || number_strategy == "Lagrange") {
    this->number_strategy_ = new NumberingLagrange< DataType >;
    this->number_strategy_->initialize(*this, this->numer_,
                                       this->numer_offsets_cell_varloc_,
                                       this->dof_interpolation_);
    this->strategy_ = true;
  } else {
    assert(0);
  }
}

/// \param description is an optional parameter that should describe what
///                    the permutation represents

template < class DataType >
void DegreeOfFreedom< DataType >::update_number_of_dofs(
  const std::string &description) {
  // Calculate number of DoFs
  this->number_of_dofs_total_ = 0;
  for (size_t i = 0, e_i = this->numer_.size(); i != e_i; ++i) {
    if (this->numer_[i] > this->number_of_dofs_total_) {
      this->number_of_dofs_total_ = this->numer_[i];
    }
  }
  ++this->number_of_dofs_total_;

  // Calculate number of Dofs for each variable
  this->number_of_dofs_for_var_.resize(nvar_, 0);
  for (size_t var = 0, e_var = this->nvar_; var != e_var; ++var) {
    int begin_offset = this->numer_offsets_cell_varloc_[var][0];
    int end_offset;

    if (var < this->nvar_ - 1) {
      end_offset = this->numer_offsets_cell_varloc_[var + 1][0];
    } else {
      end_offset = this->numer_.size();
    }

    std::vector< DofID > tmp(end_offset - begin_offset);

    for (size_t i = 0, e_i = end_offset - begin_offset; i < e_i; ++i) {
      tmp[i] = this->numer_[i + begin_offset];
    }

    std::sort(tmp.begin(), tmp.end());

    std::vector< DofID >::iterator it = std::unique(tmp.begin(), tmp.end());
    int tmp_size = it - tmp.begin();

    this->number_of_dofs_for_var_[var] = tmp_size;
  }
#if 0
  int rank = -1;
  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
  if ( ( number_of_dofs_total_old != number_of_dofs_total_ ) &&
       ( rank == 0 ) )
    std::cout << "#Dofs: " << number_of_dofs_total_
              << " " << description << std::endl;
#endif
}

template < class DataType >
void DegreeOfFreedom< DataType >::set_mesh(const mesh::Mesh *mesh) {
  this->mesh_ = mesh;
  this->mesh_flag_ = true;
}

template < class DataType >
DofID DegreeOfFreedom< DataType >::mapl2g(int var, int cell_index,
    DofID local_id) const {
  return this->numer_[this->numer_offsets_cell_varloc_[var][cell_index] +
                                                                        local_id];
}

template < class DataType >
void DegreeOfFreedom< DataType >::set_fe_manager(
  FEManager< DataType > const *manager) {
  this->fe_manager_ = manager;
  this->tdim_ = this->fe_manager_->get_dim();
  this->nvar_ = this->fe_manager_->get_nb_var();
  this->fe_manager_flag_ = true;
}

template < class DataType >
FEManager< DataType > const &
DegreeOfFreedom< DataType >::get_fe_manager() const {
  return *(this->fe_manager_);
}

/// handles DoF interpolation for hanging nodes (h- and p-refinement)

template < class DataType >
DofInterpolation &DegreeOfFreedom< DataType >::dof_interpolation() {
  return this->dof_interpolation_;
}

template < class DataType >
const DofInterpolation &DegreeOfFreedom< DataType >::dof_interpolation() const {
  return this->dof_interpolation_;
}

template < class DataType >
void DegreeOfFreedom< DataType >::get_coord_on_cell(
  int var, int cell_index, std::vector< Coord > &coords) const {
  // transformation
  CellTransformation< DataType > *transformation =
    this->get_fe_manager().get_cell_transformation(cell_index);

  // finite element type
  FEType< DataType > &fe_type = *(this->fe_manager_->get_fe_on_cell(cell_index,
                                  var));

  coords.resize(fe_type.get_nb_dof_on_cell());

  // loop over DoFs
  for (size_t i = 0, e_i = fe_type.get_nb_dof_on_cell(); i != e_i; ++i) {
    std::vector< DataType > coord_ref(get_dim());

    coord_ref[0] = fe_type.get_coord_on_cell()[i][0];
    if (get_dim() >= 2) {
      coord_ref[1] = fe_type.get_coord_on_cell()[i][1];
      if (get_dim() == 3) {
        coord_ref[2] = fe_type.get_coord_on_cell()[i][2];
      }
    }
    std::vector< DataType > coord;
    coord.resize(get_dim());

    coord[0] = transformation->x(coord_ref);
    if (get_dim() >= 2) {
      coord[1] = transformation->y(coord_ref);
      if (get_dim() == 3) {
        coord[2] = transformation->z(coord_ref);
      }
    }
    coords[i] = coord;
  }
}

template < class DataType >
void DegreeOfFreedom< DataType >::get_dofs_on_cell(
  int var, int cell_index, std::vector< DofID > &ids) const {
  // finite element type
  FEType< DataType > &fe_type = *(this->fe_manager_->get_fe_on_cell(cell_index,
                                  var));

  ids.resize(fe_type.get_nb_dof_on_cell());

  // loop over DoFs
  for (size_t i = 0, e_i = fe_type.get_nb_dof_on_cell(); i != e_i; ++i) {
    ids[i] = this->mapl2g(var, cell_index, i);
  }
}

template < class DataType >
int DegreeOfFreedom< DataType >::get_nb_dofs_on_subentity(
  int var, int tdim, int cell_index) const {
  return this->fe_manager_->get_fe_on_cell(cell_index, var)
         ->get_nb_dof_on_subentity(tdim, cell_index);
}

template < class DataType >
int DegreeOfFreedom< DataType >::get_nb_dofs_on_subentity(
  int tdim, int cell_index) const {
  int result = 0;

  for (int var = 0; var != this->fe_manager_->get_nb_var(); ++var) {
    result += this->get_nb_dofs_on_subentity(var, tdim, cell_index);
  }

  return result;
}

template < class DataType >
int DegreeOfFreedom< DataType >::get_nb_dofs_on_cell(int var,
    int cell_index) const {
  return this->fe_manager_->get_fe_on_cell(cell_index, var)->get_nb_dof_on_cell();
}

template < class DataType >
int DegreeOfFreedom< DataType >::get_nb_dofs_on_cell(int cell_index) const {
  int result = 0;

  for (int var = 0; var < this->fe_manager_->get_nb_var(); ++var) {
    result += this->get_nb_dofs_on_cell(var, cell_index);
  }

  return result;
}

template < class DataType >
void DegreeOfFreedom< DataType >::get_coord_on_subentity(
  int var, int cell_index, int tdim, int sindex,
  std::vector< Coord > &coords) const {
  // transformation
  CellTransformation< DataType > *transformation =
    this->get_fe_manager().get_cell_transformation(cell_index);

  // finite element type
  FEType< DataType > &fe_type = *(this->fe_manager_->get_fe_on_cell(cell_index,
                                  var));

  coords.resize(fe_type.get_nb_dof_on_subentity(tdim, sindex));

  // loop over DoFs
  for (size_t i = 0, e_i = fe_type.get_nb_dof_on_subentity(tdim, sindex);
       i != e_i; ++i) {
    std::vector< DataType > coord_ref(this->get_dim());

    coord_ref[0] = fe_type.get_coord_on_subentity(tdim, sindex)[i][0];
    if (this->get_dim() >= 2) {
      coord_ref[1] = fe_type.get_coord_on_subentity(tdim, sindex)[i][1];
      if (this->get_dim() == 3) {
        coord_ref[2] = fe_type.get_coord_on_subentity(tdim, sindex)[i][2];
      }
    }
    std::vector< DataType > coord;
    coord.resize(this->get_dim());

    coord[0] = transformation->x(coord_ref);
    if (this->get_dim() >= 2) {
      coord[1] = transformation->y(coord_ref);
      if (this->get_dim() == 3) {
        coord[2] = transformation->z(coord_ref);
      }
    }
    coords[i] = coord;
  }
}

template < class DataType >
void DegreeOfFreedom< DataType >::get_dofs_on_subentity(
  int var, int cell_index, int tdim, int sindex,
  std::vector< DofID > &ids) const {
  // finite element type
  FEType< DataType > &fe_type = *(this->fe_manager_->get_fe_on_cell(cell_index,
                                  var));

  ids.resize(fe_type.get_nb_dof_on_subentity(tdim, sindex));
  // loop over DoFs
  for (size_t i = 0, e_i = fe_type.get_nb_dof_on_subentity(tdim, sindex);
       i != e_i; ++i) {
    ids[i] =
      this->mapl2g(var, cell_index, fe_type.get_dof_on_subentity(tdim, sindex)[i]);
  }
}

template < class DataType >
int DegreeOfFreedom< DataType >::get_nb_dofs(int var) const {
  interminable_assert(var >= 0 && var <= this->fe_manager_->get_nb_var());
  return this->number_of_dofs_for_var_[var];
}

template < class DataType >
int DegreeOfFreedom< DataType >::get_nb_dofs() const {
  return this->number_of_dofs_total_;
}

/// \param description is an optional parameter that should describe what
///                    the permutation represents

template < class DataType >
void DegreeOfFreedom< DataType >::apply_permutation(
  const std::vector< DofID > &permutation, const std::string &description) {
  // apply permutation to mapl2g
  //
  // DoF IDs are used in numer_ only
  for (size_t i = 0, e_i = this->numer_.size(); i != e_i; ++i) {
    this->numer_[i] = permutation[this->numer_[i]];
  }

  // apply permutation to DofInterpolation
  this->dof_interpolation().apply_permutation(permutation);

  // calculate number of dofs, as this could have changed
  this->update_number_of_dofs(description);
}

/// prints numer_ to outstream

template < class DataType >
void DegreeOfFreedom< DataType >::print_numer() const {
  for (size_t i = 0, e_i = numer_.size(); i != e_i; ++i) {
    std::cout << i << "\t ->    " << numer_[i] << std::endl;
  }
}

template < class DataType >
void DegreeOfFreedom< DataType >::number(DOF_ORDERING order) {
  // check that vertex numbering of all cells is valid
  assert(this->check_mesh() == true);

  if (!this->strategy_) {
    this->init_numbering_strategy();
  }
  assert(this->strategy_);

  this->number_strategy_->number(order);
  this->update_number_of_dofs();
}

template < class DataType > int DegreeOfFreedom< DataType >::get_dim() const {
  return this->fe_manager_->get_dim();
}

template < class DataType >
void permute_constrained_dofs_to_end(DegreeOfFreedom< DataType > &dof) {
  const int num_dofs = dof.get_nb_dofs();
  std::vector< int > permutation(num_dofs, -1);
  for (int i = 0; i != num_dofs; ++i) {
    permutation[i] = i;
  }

  const DofInterpolation &dof_interp = dof.dof_interpolation();

  // Compute permutation that places constrained dofs after
  // unconstrained dofs. This is accomplished by keeping two
  // references, one to the first constrained dof, and one to the
  // last unconstrained dof. These are updated iteratively, and the
  // positions in the permutation are swapped until
  int first_constrained = 0;
  int last_unconstrained = num_dofs - 1;

  while (true) {
    // make first_constrained point to first constrained dof
    while (first_constrained < last_unconstrained &&
           dof_interp.count(first_constrained) == 0) {
      ++first_constrained;
    }

    // make last_unconstrained point to last unconstrained dof
    while (first_constrained < last_unconstrained &&
           dof_interp.count(last_unconstrained) == 1) {
      --last_unconstrained;
    }

    // if we are not done, swap the two positions
    if (first_constrained < last_unconstrained) {
      std::swap(permutation[first_constrained],
                permutation[last_unconstrained]);

      // update pointers here
      ++first_constrained;
      --last_unconstrained;
    } else {
      // we are done, break out of loop
      break;
    }
  }

  // Apply the permutation
  dof.apply_permutation(permutation);
}

template class DegreeOfFreedom< double >;
template class DegreeOfFreedom< float >;

template void permute_constrained_dofs_to_end(DegreeOfFreedom< double > &dof);
template void permute_constrained_dofs_to_end(DegreeOfFreedom< float > &dof);
} // namespace doffem
} // namespace hiflow
