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

#ifndef _DOF_DOF_PARTITION_LOCAL_H_
#define _DOF_DOF_PARTITION_LOCAL_H_

#include <map>
#include <vector>

#include "common/log.h"
#include "common/vector_algebra.h"
#include "dof/dof_interpolation.h"

namespace hiflow {
namespace mesh {
class Mesh;
}

namespace doffem {

template < class DataType, int DIM > class FEManager;
template < class DataType, int DIM > class RefElement;
template < class DataType > class DofInterpolation;

/// \author Michael Schick<br>Martin Baumann<br>Philipp Gerstner

template < class DataType, int DIM > 
class DofPartitionLocal 
{
public:
  /// Constructor
  DofPartitionLocal();

  /// Destructor
  virtual ~DofPartitionLocal();

  /// Set given mesh to a constant pointer
  void set_mesh(const mesh::Mesh *mesh);

  /// \brief Mapping local 2 global. Local DofId is local to a specific
  ///        cell and a specific variable. 
  /// Global DofId is the DOF Id on the local mesh, unless DofPartitonGlobal::renumber is called.
  /// Then, global id refers to the whole mesh
    
  DofID cell2global(size_t fe_ind, int cell_index, DofID local_id) const;
  DataType cell2factor(size_t fe_ind, int cell_index, DofID local_id) const;

  /// Setting the FEManager for the given mesh
  void set_fe_manager(FEManager< DataType, DIM > const *manager);

  /// Get the FEManager for the given mesh
  FEManager< DataType, DIM > const &get_fe_manager() const;

  /// Get the mesh
  const mesh::Mesh &get_mesh() const { return *mesh_; }

  /// Get the DoF Interpolation
  /// handles DoF interpolation for hanging nodes (h- and p-refinement)
  DofInterpolation<DataType> &dof_interpolation();
  const DofInterpolation<DataType> &dof_interpolation() const;

  /// Get the global DofIds (w.r.t. complete mesh) on a specific mesh cell and variable
  void get_dofs_on_cell(size_t fe_ind, int cell_index, std::vector< DofID > &ids) const;
  
  /// Get the number of dofs for a specific variable on a specific cell
  size_t nb_dofs_on_cell(size_t fe_ind, int cell_index) const;
  
  /// Get the total number of dofs on a specific cell (including all variables)
  size_t nb_dofs_on_cell(int cell_index) const;

  /// Get the number of dofs on the boundary for a specific variable on a
  /// specific cell
  size_t nb_dofs_on_subentity(size_t fe_ind, int tdim, int cell_index) const;
  
  /// Get the total number of dofs on the boundary on a specific cell (including
  /// all variables)
  size_t nb_dofs_on_subentity(int tdim, int cell_index) const;

  /// Get the global DofIds on a subentity (point, edge, fase) of a cell
  void get_dofs_on_subentity(size_t fe_ind, int cell_index, int tdim, int sindex, std::vector< DofID > &ids) const;
  
  /// Get the number of dofs on the local subdomain for a specific variable
  size_t nb_dofs_local(size_t fe_ind) const;
  
  /// Get the number of dofs on the subdomain for all variables
  size_t nb_dofs_local() const;

  /// Apply permutation of DoF IDs
  /// \param description is an optional parameter that should describe what
  ///                    the permutation represents
  void apply_permutation(const std::vector< DofID > &permutation, const std::string & = "");

  /// Update number_of_dofs_total_ and number_of_dofs_for_var_
  /// \param description is an optional parameter that should describe what
  ///                    the permutation represents
  void update_number_of_dofs(const std::string & = "");

  /// Printing information about the numbering field
  void print_numer() const;

  void set_applied_number_strategy (bool flag) 
  {
    this->applied_number_strategy_ = flag;
  }

  std::vector< DofID >* numer_cell_2_subdom() 
  {
    return &numer_cell_2_global_;
  }
  
  std::vector< DataType >* numer_cell_2_factor() 
  {
    return &numer_cell_2_factor_;
  }

  std::vector< std::vector< int > >* numer_cell_2_global_offsets()
  {
    return &numer_cell_2_global_offsets_;
  }

  /// Check that ordering of vertices of mesh entity is valid.
  bool check_mesh();
  
  inline std::vector<int> get_phys_var(int dof_var) const;

  inline size_t var_2_fe(size_t var) const
  {
    assert (this->fe_manager_ != NULL);
    return this->fe_manager_->var_2_fe (var);
  }

  inline size_t nb_fe () const 
  {
    return this->nb_fe_;
  }

  inline size_t nb_var () const 
  {
    return this->nb_var_;
  }

  inline size_t tdim() const 
  {
    return this->tdim_;
  }

protected:
  /// Topological dimension
  size_t tdim_;

  /// Total number of variables
  size_t nb_var_;
  size_t nb_fe_;

  /// Holds DoF IDs, needed for cell2global
  std::vector< DofID > numer_cell_2_global_;
  std::vector< DataType > numer_cell_2_factor_; 
  
  /// Offset container for numer_, needed for cell2global
  std::vector< std::vector< int > > numer_cell_2_global_offsets_;

  /// Const pointer to mesh
  const mesh::Mesh *mesh_;

  /// FEManager on the given mesh
  FEManager< DataType, DIM > const *fe_manager_;

  /// Interpolation Container, which stores the interpolation weights
  DofInterpolation<DataType> dof_interpolation_;

  /// Total number of dofs for all variables
  int local_nb_dofs_total_;

  /// Total number of dofs per variable
  std::vector< size_t > local_nb_dofs_for_fe_;

  bool applied_number_strategy_;
  
private:
  /// Check if mesh is set
  bool mesh_flag_;
  /// Check if fe_manager is set
  bool fe_manager_flag_;
};


} // namespace doffem
} // namespace hiflow
#endif
