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

#ifndef _DOF_DEGREE_OF_FREEDOM_H_
#define _DOF_DEGREE_OF_FREEDOM_H_

#include <map>
#include <vector>

#include "common/vector_algebra.h"
#include "dof/dof_fem_types.h"
#include "dof/dof_interpolation.h"
#include "fem/cell_transformation.h"
#include "fem/femanager.h"
#include "mesh/mesh.h"
#include "mesh/iterator.h"
#include "mesh/cell_type.h"

namespace hiflow {
namespace doffem {

template < class DataType, int DIM > class FEManager;

template < class DataType, int DIM > class FEType;

/// \author Michael Schick<br>Martin Baumann

template < class DataType, int DIM > class DegreeOfFreedom {
public:
  typedef Vec<DIM, DataType> Coord;

  /// Constructor
  DegreeOfFreedom();

  /// Destructor
  virtual ~DegreeOfFreedom();

  /// Set given mesh to a constant pointer
  void set_mesh(const mesh::Mesh *mesh);

  /// \brief Mapping local 2 global. Local DofId is local to a specific
  ///        cell and a specific variable. Global DofId is the global Id
  ///        on the whole given mesh
  DofID aaa_mapl2g(int fe_var, int cell_index, DofID local_id) const;

  /// Setting the FEManager for the given mesh
  void set_fe_manager(FEManager< DataType, DIM > const *manager);

  /// Get the FEManager for the given mesh
  FEManager< DataType, DIM > const &get_fe_manager() const;

  /// Get the mesh

  const mesh::Mesh &get_mesh() const { return *mesh_; }

  /// Get the DoF Interpolation
  DofInterpolation &dof_interpolation();
  const DofInterpolation &dof_interpolation() const;

  /// Get the coordinates of dofs on a mesh cell for a specific variable and
  /// cell index
  void aaa_get_dof_coord_on_cell(int dof_var, int cell_index, std::vector< Coord > &coords) const;
  void bbb_get_dof_coord_on_cell(int phys_var, int cell_index, std::vector< Coord > &coords) const 
  {
    this->aaa_get_dof_coord_on_cell(this->get_dof_var(phys_var), cell_index, coords);
  }

  void aaa_get_dof_coord_on_cell(int dof_var, int cell_index, std::vector< DataType > &coords) const;
  void bbb_get_dof_coord_on_cell(int phys_var, int cell_index, std::vector< DataType > &coords) const 
  {
    this->aaa_get_dof_coord_on_cell(this->get_dof_var(phys_var), cell_index, coords);
  }

  /// Get the global DofIds on a specific mesh cell and variable
  void aaa_get_dofs_on_cell(int dof_var, int cell_index, std::vector< DofID > &ids) const;
  void bbb_get_dofs_on_cell(int phys_var, int cell_index, std::vector< DofID > &ids) const 
  {
    this->aaa_get_dofs_on_cell(this->get_dof_var(phys_var), cell_index, ids);
  }

  /// Get the number of dofs for a specific variable on a specific cell
  int aaa_get_nb_dofs_on_cell(int dof_var, int cell_index) const;
  int bbb_get_nb_dofs_on_cell(int phys_var, int cell_index) const 
  {
    return this->aaa_get_nb_dofs_on_cell(this->get_dof_var(phys_var), cell_index);
  }

  /// Get the total number of dofs on a specific cell (including all variables)
  int get_nb_dofs_on_cell(int cell_index) const;

  /// Get the number of dofs on the boundary for a specific variable on a
  /// specific cell
  int aaa_get_nb_dofs_on_subentity(int dof_var, int tdim, int cell_index) const;
  int bbb_get_nb_dofs_on_subentity(int phys_var, int tdim, int cell_index) const 
  {
    return this->aaa_get_nb_dofs_on_subentity(this->get_dof_var(phys_var), tdim, cell_index);
  }

  /// Get the total number of dofs on the boundary on a specific cell (including
  /// all variables)
  int get_nb_dofs_on_subentity(int tdim, int cell_index) const;

  /// Get the coordinates of dofs on a subentity (point, edge, fase) of a cell
  void aaa_get_dof_coord_on_subentity(int dof_var, int cell_index, int tdim, int sindex, std::vector< Coord > &coords) const;
  void bbb_get_dof_coord_on_subentity(int phys_var, int cell_index, int tdim, int sindex, std::vector< Coord > &coords) const 
  {
    this->aaa_get_dof_coord_on_subentity(this->get_dof_var(phys_var), cell_index, tdim, sindex, coords);
  }

  /// Get the global DofIds on a subentity (point, edge, fase) of a cell
  void aaa_get_dofs_on_subentity(int dof_var, int cell_index, int tdim, int sindex, std::vector< DofID > &ids) const;
  void bbb_get_dofs_on_subentity(int phys_var, int cell_index, int tdim, int sindex, std::vector< DofID > &ids) const 
  {
    this->aaa_get_dofs_on_subentity(this->get_dof_var(phys_var), cell_index, tdim, sindex, ids);
  }

  /// Get the number of dofs on the whole mesh for a specific variable
  int aaa_get_nb_dofs(int dof_var) const;
  int bbb_get_nb_dofs(int phys_var) const 
  {
    return this->aaa_get_nb_dofs(this->get_dof_var(phys_var));
  }

  /// Get the number of dofs on the whole mesh for all variables
  int get_nb_dofs() const;

  /// Apply permutation of DoF IDs
  void apply_permutation(const std::vector< DofID > &permutation, const std::string & = "");

  /// Update number_of_dofs_total_ and number_of_dofs_for_var_
  void update_number_of_dofs(const std::string & = "");

  /// Printing information about the numbering field
  void print_numer() const;

  void set_applied_number_strategy (bool flag) 
  {
    this->applied_number_strategy_ = flag;
  }

  std::vector< DofID >* numer() 
  {
    return &numer_;
  }

  std::vector< std::vector< int > >* numer_offsets_cell_dofvarloc()
  {
    return &numer_offsets_cell_dofvarloc_;
  }

  /// Check that ordering of vertices of mesh entity is valid.
  bool check_mesh();
  
  inline std::vector<int> get_phys_var(int dof_var) const;

  inline int get_dof_var(int phys_var) const;

protected:
  /// Topological dimension
  int tdim_;
  /// Total number of variables
  int dof_nvar_;

  /// Holds DoF IDs, needed for mapl2g
  std::vector< DofID > numer_;

  /// Offset container for numer_, needed for mapl2g
  std::vector< std::vector< int > > numer_offsets_cell_dofvarloc_;

  /// Const pointer to mesh
  const mesh::Mesh *mesh_;

  /// FEManager on the given mesh
  FEManager< DataType, DIM > const *fe_manager_;

  /// Interpolation Container, which stores the interpolation weights
  DofInterpolation dof_interpolation_;

  /// Get the geometrical dimension stored in FEManager
  int get_dim() const;

  /// Total number of dofs for all variables
  int number_of_dofs_total_;

  /// Total number of dofs per variable
  std::vector< int > number_of_dofs_for_var_;

  bool applied_number_strategy_;
  
private:
  /// Check if mesh is set
  bool mesh_flag_;
  /// Check if fe_manager is set
  bool fe_manager_flag_;
};

template < class DataType, int DIM > 
bool DegreeOfFreedom< DataType, DIM >::check_mesh() {
  bool ret = true;

  // loop over mesh cells
  for (mesh::EntityIterator it = mesh_->begin(get_dim()),
                            e_it = mesh_->end(get_dim());
       it != e_it; ++it) {
    // get coordinates of current cell
    std::vector< mesh::Coordinate > coords;
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

template < class DataType, int DIM > 
DegreeOfFreedom< DataType, DIM >::DegreeOfFreedom() {
  mesh_ = NULL;

  number_of_dofs_total_ = -1;

  mesh_flag_ = false;
  fe_manager_flag_ = false;

  this->applied_number_strategy_ = false;
}

template < class DataType, int DIM > 
DegreeOfFreedom< DataType, DIM >::~DegreeOfFreedom() {
}

/// \param description is an optional parameter that should describe what
///                    the permutation represents

template < class DataType, int DIM >
void DegreeOfFreedom< DataType, DIM >::update_number_of_dofs(
    const std::string &description) {
  // Calculate number of DoFs
  number_of_dofs_total_ = 0;
  for (size_t i = 0, e_i = numer_.size(); i != e_i; ++i) {
    if (numer_[i] > number_of_dofs_total_) {
      number_of_dofs_total_ = numer_[i];
    }
  }
  ++number_of_dofs_total_;

  // Calculate number of Dofs for each variable
  number_of_dofs_for_var_.resize(dof_nvar_, 0);
  for (size_t var = 0, e_var = dof_nvar_; var != e_var; ++var) {
    int begin_offset = numer_offsets_cell_dofvarloc_[var][0];
    int end_offset;

    if (var < dof_nvar_ - 1) {
      end_offset = numer_offsets_cell_dofvarloc_[var + 1][0];
    } else {
      end_offset = numer_.size();
    }

    std::vector< DofID > tmp(end_offset - begin_offset);

    for (size_t i = 0, e_i = end_offset - begin_offset; i < e_i; ++i) {
      tmp[i] = numer_[i + begin_offset];
    }

    std::sort(tmp.begin(), tmp.end());

    std::vector< DofID >::iterator it = std::unique(tmp.begin(), tmp.end());
    int tmp_size = it - tmp.begin();

    number_of_dofs_for_var_[var] = tmp_size;
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

template < class DataType, int DIM >
void DegreeOfFreedom< DataType, DIM >::set_mesh(const mesh::Mesh *mesh) {
  mesh_ = mesh;
  mesh_flag_ = true;
}

template < class DataType, int DIM >
DofID DegreeOfFreedom< DataType, DIM >::aaa_mapl2g(int dof_var, int cell_index,
                                          DofID local_id) const {
  return numer_[numer_offsets_cell_dofvarloc_[dof_var][cell_index] + local_id];
}

template < class DataType, int DIM >
void DegreeOfFreedom< DataType, DIM >::set_fe_manager(
    FEManager< DataType, DIM > const *manager) {
  fe_manager_ = manager;
  tdim_ = fe_manager_->get_dim();
  dof_nvar_ = fe_manager_->aaa_get_nb_var();
  fe_manager_flag_ = true;
}

template < class DataType, int DIM >
FEManager< DataType, DIM > const & DegreeOfFreedom< DataType, DIM >::get_fe_manager() const {
  return *fe_manager_;
}

/// handles DoF interpolation for hanging nodes (h- and p-refinement)

template < class DataType, int DIM >
DofInterpolation &DegreeOfFreedom< DataType, DIM >::dof_interpolation() {
  return dof_interpolation_;
}

template < class DataType, int DIM >
const DofInterpolation &DegreeOfFreedom< DataType, DIM >::dof_interpolation() const {
  return dof_interpolation_;
}

template < class DataType, int DIM >
void DegreeOfFreedom< DataType, DIM >::aaa_get_dof_coord_on_cell( int dof_var, int cell_index, std::vector< Coord > &coords) const {
  // transformation
  CellTransformation< DataType, DIM > *transformation = get_fe_manager().get_cell_transformation(cell_index);

  // finite element type
  FEType< DataType, DIM > &fe_type = *(fe_manager_->aaa_get_fe_on_cell(cell_index, dof_var));

  coords.resize(fe_type.get_nb_dof_on_cell());

  // loop over DoFs
  for (size_t i = 0, e_i = fe_type.get_nb_dof_on_cell(); i != e_i; ++i) {

    coords[i][0] = transformation->x(fe_type.get_dof_coord_on_cell()[i]);
    if (get_dim() >= 2) {
      coords[i][1] = transformation->y(fe_type.get_dof_coord_on_cell()[i]);
      if (get_dim() == 3) {
        coords[i][2] = transformation->z(fe_type.get_dof_coord_on_cell()[i]);
      }
    }
  }
}

template < class DataType, int DIM >
void DegreeOfFreedom< DataType, DIM >::aaa_get_dof_coord_on_cell( int dof_var, int cell_index, std::vector< DataType > &coords) const {
  // transformation
  CellTransformation< DataType, DIM > *transformation = get_fe_manager().get_cell_transformation(cell_index);

  // finite element type
  FEType< DataType, DIM > &fe_type = *(fe_manager_->aaa_get_fe_on_cell(cell_index, dof_var));

  coords.resize(fe_type.get_nb_dof_on_cell() * DIM);

  // loop over DoFs
  for (size_t i = 0, e_i = fe_type.get_nb_dof_on_cell(); i != e_i; ++i) {

    coords[i * DIM + 0] = transformation->x(fe_type.get_dof_coord_on_cell()[i]);
    if (get_dim() >= 2) {
      coords[i * DIM + 1] = transformation->y(fe_type.get_dof_coord_on_cell()[i]);
      if (get_dim() == 3) {
        coords[i * DIM + 2] = transformation->z(fe_type.get_dof_coord_on_cell()[i]);
      }
    }
  }
}

template < class DataType, int DIM >
void DegreeOfFreedom< DataType, DIM >::aaa_get_dofs_on_cell(
    int dof_var, int cell_index, std::vector< DofID > &ids) const {
  // finite element type
  FEType< DataType, DIM > &fe_type = *(fe_manager_->aaa_get_fe_on_cell(cell_index, dof_var));

  ids.resize(fe_type.get_nb_dof_on_cell());

  // loop over DoFs
  for (size_t i = 0, e_i = fe_type.get_nb_dof_on_cell(); i != e_i; ++i) {
    ids[i] = aaa_mapl2g(dof_var, cell_index, i);
  }
}

template < class DataType, int DIM >
int DegreeOfFreedom< DataType, DIM >::aaa_get_nb_dofs_on_subentity(
    int dof_var, int tdim, int cell_index) const {
  return fe_manager_->aaa_get_fe_on_cell(cell_index, dof_var)
      ->get_nb_dof_on_subentity(tdim, cell_index);
}

template < class DataType, int DIM >
int DegreeOfFreedom< DataType, DIM >::get_nb_dofs_on_subentity(
    int tdim, int cell_index) const {
  int result = 0;

  for (int var = 0; var < fe_manager_->aaa_get_nb_var(); ++var) {
    result += aaa_get_nb_dofs_on_subentity(var, tdim, cell_index);
  }

  return result;
}

template < class DataType, int DIM >
int DegreeOfFreedom< DataType, DIM >::aaa_get_nb_dofs_on_cell(int dof_var,
                                                     int cell_index) const {
  return fe_manager_->aaa_get_fe_on_cell(cell_index, dof_var)->get_nb_dof_on_cell();
}

template < class DataType, int DIM >
int DegreeOfFreedom< DataType, DIM >::get_nb_dofs_on_cell(int cell_index) const {
  int result = 0;

  for (int var = 0; var < fe_manager_->aaa_get_nb_var(); ++var) {
    result += aaa_get_nb_dofs_on_cell(var, cell_index);
  }

  return result;
}

template < class DataType, int DIM >
void DegreeOfFreedom< DataType, DIM >::aaa_get_dof_coord_on_subentity(
    int dof_var, int cell_index, int tdim, int sindex,
    std::vector< Coord > &coords) const {
  // transformation
  CellTransformation< DataType, DIM > *transformation =
      get_fe_manager().get_cell_transformation(cell_index);

  // finite element type
  FEType< DataType, DIM > &fe_type = *(fe_manager_->aaa_get_fe_on_cell(cell_index, dof_var));

  coords.resize(fe_type.get_nb_dof_on_subentity(tdim, sindex));
  // loop over DoFs
  for (size_t i = 0, e_i = fe_type.get_nb_dof_on_subentity(tdim, sindex); i != e_i; ++i) {

    coords[i][0] = transformation->x(fe_type.get_dof_coord_on_subentity(tdim, sindex)[i]);
    if (get_dim() >= 2) {
      coords[i][1] = transformation->y(fe_type.get_dof_coord_on_subentity(tdim, sindex)[i]);
      if (get_dim() == 3) {
        coords[i][2] = transformation->z(fe_type.get_dof_coord_on_subentity(tdim, sindex)[i]);
      }
    }
  }
}

template < class DataType, int DIM >
void DegreeOfFreedom< DataType, DIM >::aaa_get_dofs_on_subentity(
    int dof_var, int cell_index, int tdim, int sindex,
    std::vector< DofID > &ids) const {
  // finite element type
  FEType< DataType, DIM > &fe_type = *(fe_manager_->aaa_get_fe_on_cell(cell_index, dof_var));

  ids.resize(fe_type.get_nb_dof_on_subentity(tdim, sindex));
  // loop over DoFs
  for (size_t i = 0, e_i = fe_type.get_nb_dof_on_subentity(tdim, sindex);
       i != e_i; ++i) {
    ids[i] =
        aaa_mapl2g(dof_var, cell_index, fe_type.get_dof_on_subentity(tdim, sindex)[i]);
  }
}

template < class DataType, int DIM >
int DegreeOfFreedom< DataType, DIM >::aaa_get_nb_dofs(int dof_var) const {
  interminable_assert(dof_var >= 0 && dof_var <= fe_manager_->aaa_get_nb_var());
  return number_of_dofs_for_var_[dof_var];
}

template < class DataType, int DIM >
int DegreeOfFreedom< DataType, DIM >::get_nb_dofs() const {
  return number_of_dofs_total_;
}

template < class DataType, int DIM >
int DegreeOfFreedom< DataType, DIM >::get_dof_var(int phys_var) const {
  assert (this->fe_manager_ != NULL);
  return this->fe_manager_->get_dof_var (phys_var);
}

template < class DataType, int DIM >
std::vector<int> DegreeOfFreedom< DataType, DIM >::get_phys_var(int dof_var) const {
  assert (this->fe_manager_ != NULL);
  return this->fe_manager_->get_phys_var (dof_var);
}

/// \param description is an optional parameter that should describe what
///                    the permutation represents

template < class DataType, int DIM >
void DegreeOfFreedom< DataType, DIM >::apply_permutation(
    const std::vector< DofID > &permutation, const std::string &description) {
  // apply permutation to mapl2g
  //
  // DoF IDs are used in numer_ only
  for (size_t i = 0, e_i = numer_.size(); i != e_i; ++i) {
    numer_[i] = permutation[numer_[i]];
  }

  // apply permutation to DofInterpolation
  dof_interpolation().apply_permutation(permutation);

  // calculate number of dofs, as this could have changed
  update_number_of_dofs(description);
}

/// prints numer_ to outstream

template < class DataType, int DIM >
void DegreeOfFreedom< DataType, DIM >::print_numer() const {
  for (size_t i = 0, e_i = numer_.size(); i != e_i; ++i) {
    std::cout << i << "\t ->    " << numer_[i] << std::endl;
  }
}

template < class DataType, int DIM > int DegreeOfFreedom< DataType, DIM >::get_dim() const {
  return fe_manager_->get_dim();
}

template < class DataType, int DIM >
void permute_constrained_dofs_to_end(DegreeOfFreedom< DataType, DIM > &dof) {
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


} // namespace doffem
} // namespace hiflow
#endif
