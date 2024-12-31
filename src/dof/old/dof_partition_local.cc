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

#include "dof_partition_local.h"
#include "dof/dof_interpolation.h"
#include "fem/fe_manager.h"
#include "fem/fe_reference.h"
#include "mesh/mesh.h"
#include "mesh/iterator.h"
#include "mesh/cell_type.h"

namespace hiflow {
namespace doffem {

template < class DataType, int DIM > 
bool DofPartitionLocal< DataType, DIM >::check_mesh() 
{
  bool ret = true;

  // loop over mesh cells
  for (mesh::EntityIterator it = mesh_->begin(tdim()), e_it = mesh_->end(tdim()); it != e_it; ++it) 
  {
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
DofPartitionLocal< DataType, DIM >::DofPartitionLocal() 
{
  mesh_ = NULL;

  local_nb_dofs_total_ = -1;

  mesh_flag_ = false;
  fe_manager_flag_ = false;

  this->applied_number_strategy_ = false;
}

template < class DataType, int DIM > 
DofPartitionLocal< DataType, DIM >::~DofPartitionLocal() 
{
}

template < class DataType, int DIM >
void DofPartitionLocal< DataType, DIM >::update_number_of_dofs( const std::string &description) 
{
  // Calculate number of DoFs
  local_nb_dofs_total_ = 0;
  for (size_t i = 0, e_i = numer_cell_2_global_.size(); i != e_i; ++i) 
  {
    if (numer_cell_2_global_[i] > local_nb_dofs_total_) 
    {
      local_nb_dofs_total_ = numer_cell_2_global_[i];
    }
  }
  ++local_nb_dofs_total_;

  // Calculate number of Dofs for each variable
  local_nb_dofs_for_fe_.resize(nb_fe_, 0);
  for (size_t fe_ind = 0; fe_ind != nb_fe_; ++fe_ind) 
  {
    int begin_offset = numer_cell_2_global_offsets_[fe_ind][0];
    int end_offset;

    if (fe_ind < nb_fe_ - 1) 
    {
      end_offset = numer_cell_2_global_offsets_[fe_ind + 1][0];
    } 
    else 
    {
      end_offset = numer_cell_2_global_.size();
    }

    std::vector< DofID > tmp(end_offset - begin_offset);

    for (size_t i = 0, e_i = end_offset - begin_offset; i < e_i; ++i) 
    {
      tmp[i] = numer_cell_2_global_[i + begin_offset];
    }

    std::sort(tmp.begin(), tmp.end());

    std::vector< DofID >::iterator it = std::unique(tmp.begin(), tmp.end());
    int tmp_size = it - tmp.begin();

    local_nb_dofs_for_fe_[fe_ind] = tmp_size;
  }
#if 0
            int rank = -1;
            MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
            if ( ( number_of_dofs_total_old != local_nb_dofs_total_ ) &&
                 ( rank == 0 ) )
                std::cout << "#Dofs: " << local_nb_dofs_total_
                    << " " << description << std::endl;
#endif
}

template < class DataType, int DIM >
void DofPartitionLocal< DataType, DIM >::set_mesh(const mesh::Mesh *mesh) 
{
  mesh_ = mesh;
  mesh_flag_ = true;
}

template < class DataType, int DIM >
DofID DofPartitionLocal< DataType, DIM >::cell2global(size_t fe_ind, int cell_index, DofID local_id) const 
{
  return numer_cell_2_global_[numer_cell_2_global_offsets_[fe_ind][cell_index] + local_id];
}

template < class DataType, int DIM >
DataType DofPartitionLocal< DataType, DIM >::cell2factor(size_t fe_ind, int cell_index, DofID local_id) const 
{
  return numer_cell_2_factor_[numer_cell_2_global_offsets_[fe_ind][cell_index] + local_id];
}

template < class DataType, int DIM >
void DofPartitionLocal< DataType, DIM >::set_fe_manager( FEManager< DataType, DIM > const *manager) 
{
  fe_manager_ = manager;
  tdim_ = fe_manager_->tdim();
  nb_fe_ = fe_manager_->nb_fe();
  nb_var_ = fe_manager_->nb_var();
  fe_manager_flag_ = true;
}

template < class DataType, int DIM >
FEManager< DataType, DIM > const & DofPartitionLocal< DataType, DIM >::get_fe_manager() const 
{
  return *fe_manager_;
}

template < class DataType, int DIM >
DofInterpolation<DataType> &DofPartitionLocal< DataType, DIM >::dof_interpolation() 
{
  return dof_interpolation_;
}

template < class DataType, int DIM >
const DofInterpolation<DataType> &DofPartitionLocal< DataType, DIM >::dof_interpolation() const 
{
  return dof_interpolation_;
}

template < class DataType, int DIM >
void DofPartitionLocal< DataType, DIM >::get_dofs_on_cell( size_t fe_ind, int cell_index, std::vector< DofID > &ids) const 
{
  // finite element type
  const RefElement< DataType, DIM > &ref_fe = *(fe_manager_->get_fe(cell_index, fe_ind));

  ids.resize(ref_fe.nb_dof_on_cell());

  // loop over DoFs
  for (size_t i = 0, e_i = ref_fe.nb_dof_on_cell(); i != e_i; ++i) 
  {
    ids[i] = this->cell2global(fe_ind, cell_index, i);
  }
}

template < class DataType, int DIM >
size_t DofPartitionLocal< DataType, DIM >::nb_dofs_on_subentity( size_t fe_ind , int tdim, int cell_index) const 
{
  return fe_manager_->get_fe(cell_index, fe_ind)->nb_dof_on_subentity(tdim, cell_index);
}

template < class DataType, int DIM >
size_t DofPartitionLocal< DataType, DIM >::nb_dofs_on_subentity( int tdim, int cell_index) const 
{
  int result = 0;

  for (size_t fe_ind = 0; fe_ind < fe_manager_->nb_fe(); ++fe_ind) 
  {
    result += this->nb_dofs_on_subentity(fe_ind, tdim, cell_index);
  }

  return result;
}

template < class DataType, int DIM >
size_t DofPartitionLocal< DataType, DIM >::nb_dofs_on_cell(size_t fe_ind, int cell_index) const 
{
  return fe_manager_->get_fe(cell_index, fe_ind)->nb_dof_on_cell();
}

template < class DataType, int DIM >
size_t DofPartitionLocal< DataType, DIM >::nb_dofs_on_cell(int cell_index) const 
{
  int result = 0;

  for (size_t fe_ind = 0; fe_ind < fe_manager_->nb_fe(); ++fe_ind)  
  {
    result += this->nb_dofs_on_cell(fe_ind, cell_index);
  }

  return result;
}

template < class DataType, int DIM >
void DofPartitionLocal< DataType, DIM >::get_dofs_on_subentity( size_t fe_ind, int cell_index, int tdim, int sindex,
                                                                std::vector< DofID > &ids) const 
{
  // finite element type
  const RefElement< DataType, DIM > &fe_type = *(fe_manager_->get_fe(cell_index, fe_ind));

  ids.resize(fe_type.nb_dof_on_subentity(tdim, sindex));
  
  // loop over DoFs
  for (size_t i = 0, e_i = fe_type.nb_dof_on_subentity(tdim, sindex); i != e_i; ++i) 
  {
    ids[i] = this->cell2global(fe_ind, cell_index, fe_type.get_dof_on_subentity(tdim, sindex)[i]);
  }
}

template < class DataType, int DIM >
size_t DofPartitionLocal< DataType, DIM >::nb_dofs_local(size_t fe_ind) const 
{
  interminable_assert(fe_ind >= 0 && fe_ind <= fe_manager_->nb_fe());
  return local_nb_dofs_for_fe_[fe_ind];
}

template < class DataType, int DIM >
size_t DofPartitionLocal< DataType, DIM >::nb_dofs_local() const 
{
  return local_nb_dofs_total_;
}

template < class DataType, int DIM >
void DofPartitionLocal< DataType, DIM >::apply_permutation( const std::vector< DofID > &permutation, 
                                                            const std::string &description) 
{
  // apply permutation to cell2global
  //
  // DoF IDs are used in numer_ only
  for (size_t i = 0, e_i = numer_cell_2_global_.size(); i != e_i; ++i) 
  {
    numer_cell_2_global_[i] = permutation[numer_cell_2_global_[i]];
  }
#if 0
CHANGE
permute numer_cell_2_factor
#endif

  // apply permutation to DofInterpolation
  dof_interpolation().apply_permutation(permutation);

  // calculate number of dofs, as this could have changed
  update_number_of_dofs(description);
}

template < class DataType, int DIM >
void DofPartitionLocal< DataType, DIM >::print_numer() const 
{ 
  for (size_t i = 0, e_i = numer_cell_2_global_.size(); i != e_i; ++i) 
  {
    std::cout << i << "\t ->    " << numer_cell_2_global_[i] << " : " << numer_cell_2_factor_[i] << std::endl;
  }
}

template < class DataType, int DIM >
void permute_constrained_dofs_to_end(DofPartitionLocal< DataType, DIM > &dof) 
{
  const int num_dofs = dof.nb_dofs();
  std::vector< int > permutation(num_dofs, -1);
  for (int i = 0; i != num_dofs; ++i) 
  {
    permutation[i] = i;
  }

  const DofInterpolation<DataType> &dof_interp = dof.dof_interpolation();

  // Compute permutation that places constrained dofs after
  // unconstrained dofs. This is accomplished by keeping two
  // references, one to the first constrained dof, and one to the
  // last unconstrained dof. These are updated iteratively, and the
  // positions in the permutation are swapped until
  int first_constrained = 0;
  int last_unconstrained = num_dofs - 1;

  while (true) 
  {
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

template class DofPartitionLocal< double, 3 >;
template class DofPartitionLocal< double, 2 >;
template class DofPartitionLocal< double, 1 >;

template class DofPartitionLocal< float, 3 >;
template class DofPartitionLocal< float, 2 >;
template class DofPartitionLocal< float, 1 >;

} // namespace doffem
} // namespace hiflow
