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

#include "dof/dof_partition_global.h"
#include "mesh/mesh.h"
#include "mesh/iterator.h"
#include "fem/fe_manager.h"
#include "fem/fe_reference.h"

namespace hiflow {
namespace doffem {

template < class DataType, int DIM >
const MPI_Comm &DofPartitionGlobal< DataType, DIM >::get_mpi_comm() const 
{
  return comm_;
}

template < class DataType, int DIM >
int DofPartitionGlobal< DataType, DIM >::my_subdom() const 
{
  return my_subdom_;
}

template < class DataType, int DIM >
int DofPartitionGlobal< DataType, DIM >::my_dof_offset() const 
{
  return my_dof_offset_;
}

template < class DataType, int DIM >
size_t DofPartitionGlobal< DataType, DIM >::nb_subdom() const 
{
  return this->nb_subdom_;
}

template < class DataType, int DIM >
size_t DofPartitionGlobal< DataType, DIM >::nb_dofs_on_subdom(int subdomain) const 
{
  assert(subdomain >= 0 && subdomain < nb_subdom_);
  return nb_dofs_on_subdom_[subdomain];
}

template < class DataType, int DIM > 
size_t DofPartitionGlobal< DataType, DIM >::nb_dofs_global() const 
{
  return gl_nb_dofs_total_;
}

template < class DataType, int DIM >
void DofPartitionGlobal< DataType, DIM >::local2global(DofID local_id, DofID *global_id) const 
{
  assert(local_id >= 0 && local_id < nb_dofs_on_subdom_[my_subdom_]);
  *global_id = local_id + my_dof_offset_;
  //  *global_id = local2global_[local_id];
}

template < class DataType, int DIM >
void DofPartitionGlobal< DataType, DIM >::global2local(DofID global_id,
                                            DofID *local_id) const {
  assert(global_id >= my_dof_offset_);
  assert(global_id < my_dof_offset_ + nb_dofs_on_subdom_[my_subdom_]);
  *local_id = global_id - my_dof_offset_;
  //   std::map<DofID,DofID>::const_iterator it = global2local_.find(global_id);
  //   *local_id = (*it).second;
}

template < class DataType, int DIM >
bool DofPartitionGlobal< DataType, DIM >::is_dof_on_subdom(DofID global_id) const {
  return global_id >= my_dof_offset_ &&
         global_id < my_dof_offset_ + nb_dofs_on_subdom_[my_subdom_];
}

template < class DataType, int DIM >
int DofPartitionGlobal< DataType, DIM >::owner_of_dof(DofID global_id) const {
  int result = 0;
  int bound = 0;

  for (size_t s = 0, e_s = nb_subdom_; s != e_s; ++s) {
    bound += nb_dofs_on_subdom_[s];
    if (global_id < bound) {
      result = static_cast< int >(s);
      break;
    }
  }

  return result;
}

template < class DataType, int DIM >
DofPartitionGlobal< DataType, DIM >::DofPartitionGlobal()
    : DofPartitionLocal< DataType, DIM >::DofPartitionLocal() {
  my_dof_offset_ = 0;

  // Default Value is COMM WORLD
  comm_ = MPI_COMM_WORLD;
  MPI_Comm_size(comm_, &nb_subdom_);
  nb_dofs_on_subdom_.resize(nb_subdom_, -1);
  MPI_Comm_rank(comm_, &my_subdom_);
  shared_subdomains_.resize(nb_subdom_, false);
}

template < class DataType, int DIM > DofPartitionGlobal< DataType, DIM >::~DofPartitionGlobal() {}

template < class DataType, int DIM >
void DofPartitionGlobal< DataType, DIM >::set_mpi_comm(const MPI_Comm &comm) {
  comm_ = comm;
  MPI_Comm_size(comm, &nb_subdom_);

  nb_dofs_on_subdom_.resize(nb_subdom_, -1);
  MPI_Comm_rank(comm, &my_subdom_);
  shared_subdomains_.resize(nb_subdom_, false);
}

template < class DataType, int DIM >
void DofPartitionGlobal< DataType, DIM >::number_parallel() 
{
  assert (this->applied_number_strategy_);
  
  // Check if sequential case or truely parallel case is used
  if (nb_subdom_ == 1) {
    nb_dofs_on_subdom_[0] = this->nb_dofs_local();
  } else {
    create_ownerships();
    renumber();

    this->local_nb_dofs_total_ = nb_dofs_on_subdom_[my_subdom_];
  }

  // Calculate number of dofs on the global domain
  gl_nb_dofs_total_ = 0;
  for (size_t s = 0, e_s = nb_subdom_; s != e_s; ++s) {
    gl_nb_dofs_total_ += nb_dofs_on_subdom_[s];
  }

  // Some essential outputs
  if (my_subdom_ == 0) {
    int min_num_dof = 1e6;
    int max_num_dof = 0;
    for (size_t s = 0, e_s = nb_subdom_; s != e_s; ++s) {
      //                    LOG_INFO ( "#Dofs on subdomain " << s,
      //                    nb_dofs_on_subdom_[s] );
      if (nb_dofs_on_subdom_[s] > max_num_dof) {
        max_num_dof = nb_dofs_on_subdom_[s];
      }
      if (nb_dofs_on_subdom_[s] < min_num_dof) {
        min_num_dof = nb_dofs_on_subdom_[s];
      }
    }
    LOG_INFO("#Dofs on global domain", gl_nb_dofs_total_);
    LOG_INFO("Load balancing", static_cast< DataType >(max_num_dof) /
                                   static_cast< DataType >(min_num_dof));
  }
}

template < class DataType, int DIM > 
void DofPartitionGlobal< DataType, DIM >::create_ownerships() 
{
  ownership_.clear();
  ownership_.resize(this->nb_dofs_local(), -1); 

  for (mesh::EntityIterator it = this->mesh_->begin(this->tdim_),
                            e_it = this->mesh_->end(this->tdim_); it != e_it; ++it) 
  {
    int subdomain_of_entity;
    it->get("_sub_domain_", &subdomain_of_entity);
  
#ifndef NDEBUG
    if (subdomain_of_entity != my_subdom_) 
    {
      LOG_DEBUG(3,"[" << my_subdom_ << "] remote cell with id " << it->id());
    }
    else
    {
      LOG_DEBUG(3, "[" << my_subdom_ << "] local cell with id " << it->id());
    }
#endif
    
    for (size_t fe_ind = 0; fe_ind != this->nb_fe_; ++fe_ind) 
    {
      std::vector< DofID > local_dofs;
      this->get_dofs_on_cell(fe_ind, it->index(), local_dofs);

      for (size_t i = 0, e_i = local_dofs.size(); i != e_i; ++i) 
      {
        const DofID ld_i = local_dofs[i];
        if (subdomain_of_entity < ownership_[ld_i] || ownership_[ld_i] == -1) 
        {
          ownership_[ld_i] = subdomain_of_entity;
        }
      }
    }
  }
}

template < class DataType, int DIM > 
void DofPartitionGlobal< DataType, DIM >::renumber() 
{
  // Communicate number of dofs including ghost layer
  std::vector< int > ndofs_with_ghost( static_cast< size_t >(nb_subdom_));
  ndofs_with_ghost[static_cast< size_t >(my_subdom_)] = this->nb_dofs_local();

  int ndofs_with_ghost_sent = this->nb_dofs_local();

  // Store information about number of dofs including ghost layer
  nb_dofs_incl_ghost_ = this->nb_dofs_local();

  MPI_Allgather(&ndofs_with_ghost_sent, 1, MPI_INT, vec2ptr(ndofs_with_ghost),
                1, MPI_INT, comm_);

  LOG_DEBUG(2, "[" << my_subdom_ << "]: Number of DOFs on each subdomain: "
                   << string_from_range(ndofs_with_ghost.begin(),
                                        ndofs_with_ghost.end()));
  LOG_DEBUG(3, "[" << my_subdom_ << "]: Ownership: "
                   << string_from_range(ownership_.begin(), ownership_.end()));

  // Calculate temporary dof offset
  int tmp_dof_offset = 0;
  for (size_t s = 0, e_s = static_cast< size_t >(my_subdom_); s != e_s;
       ++s) {
    tmp_dof_offset += ndofs_with_ghost[s];
  }

  // Fill first permutation to create a local consecutive numbering w.r.t.
  // tmp_dof_offset
  std::vector< int > permutation(ownership_.size());
  int dof_number = 0;

  for (size_t i = 0, e_i = ownership_.size(); i != e_i; ++i) {
    if (ownership_[i] == my_subdom_) {
      permutation[i] = dof_number + tmp_dof_offset;
      ++dof_number;
    }
  }

  LOG_DEBUG(2,
            "[" << my_subdom_
                << "]: Number of (interior) DOFs on subdomain: " << dof_number);

  LOG_DEBUG(3,
            "[" << my_subdom_ << "]: Permutation size " << permutation.size()
                << ", content: "
                << string_from_range(permutation.begin(), permutation.end()));

  int ghost_number = 0;
  for (size_t i = 0, e_i = ownership_.size(); i != e_i; ++i) {
    if (ownership_[i] != my_subdom_) {
      permutation[i] = dof_number + tmp_dof_offset + ghost_number;
      ++ghost_number;
    }
  }

  LOG_DEBUG(2,
            "[" << my_subdom_
                << "]: Number of (ghost) DOFs on subdomain: " << ghost_number);

  LOG_DEBUG(3,
            "[" << my_subdom_ << "]: Permutation size " << permutation.size()
                << ", content: "
                << string_from_range(permutation.begin(), permutation.end()));

  this->apply_permutation(permutation);

  LOG_DEBUG(2, " First permutation done ");

  // Calculate number of dofs which belong to my subdomain
  nb_dofs_on_subdom_[static_cast< size_t >(my_subdom_)] = dof_number;

  // Gather number of dofs of all subdomains on all processes
  MPI_Allgather(&dof_number, 1, MPI_INT, vec2ptr(nb_dofs_on_subdom_), 1, MPI_INT,
                comm_);
  LOG_DEBUG(2,
            "[" << my_subdom_ << "]: ndofs_on_sd "
                << string_from_range(nb_dofs_on_subdom_.begin(), nb_dofs_on_subdom_.end()));

  // Setting up data structure for communication and management of dofs
  // which belong to ghost cells on current process.
  // The data structure maps subdomain indices to maps of (ghost) cell
  // indices. For each (ghost) cell, a map of size dof_nvar_ is hold which
  // itself contains the (global) dof numbers of this cell
  std::map< int, std::map< int, std::map< int, std::vector< int > > > >
      numer_ghost;

  // Create subdomain/ghost cell structure of numer_ghost
  for (mesh::EntityIterator it = this->mesh_->begin(this->tdim_),
                            e_it = this->mesh_->end(this->tdim_);
       it != e_it; ++it) {
    int subdomain_index;
    it->get("_sub_domain_", &subdomain_index);

    if (subdomain_index != my_subdom_) {
      // Access/create (on first access) map for remote subdomain
      numer_ghost[subdomain_index];
      int ghost_cell_index;
      it->get("_remote_index_", &ghost_cell_index);
      // Access/create (on first access) map entry for ghost cell
      numer_ghost[subdomain_index][ghost_cell_index];
    }
  }

  {
    // Number of ghost cells which I share with each other process
    std::vector< int > num_ghost_cells(
        static_cast< size_t >(this->nb_subdom_), 0);
    int total_num_ghost_cells = 0;
    for (int i = 0; i < this->nb_subdom_; ++i) {
      num_ghost_cells[i] = numer_ghost[i].size();

      LOG_DEBUG(2, "[" << my_subdom_ << "]: common ghost cells with process "
                       << i << ": " << num_ghost_cells[i]);

      total_num_ghost_cells += numer_ghost[i].size();
    }

    // Number of ghost cells which other processes share with me
    std::vector< int > num_ghost_cells_others(
        static_cast< size_t >(this->nb_subdom_), 0);

    // Exchange number of ghost cells
    MPI_Alltoall(vec2ptr(num_ghost_cells), 1, MPI_INT,
                 vec2ptr(num_ghost_cells_others), 1, MPI_INT, this->comm_);

#ifndef NDEBUG
    for (size_t i = 0; i < static_cast< size_t >(this->nb_subdom_);
         ++i) {
      LOG_DEBUG(2, "[" << my_subdom_ << "]: process " << i
                       << " requests number of ghost cells: "
                       << num_ghost_cells_others[i]);
    }
#endif

    // Ghost cell indices which I need of others
    std::vector< int > ghost_indices;
    ghost_indices.reserve(static_cast< size_t >(total_num_ghost_cells));
    for (int i = 0; i < this->nb_subdom_; ++i) {
      for (typename std::map<
               int, std::map< int, std::vector< int > > >::const_iterator
               it = numer_ghost[i].begin(),
               e_it = numer_ghost[i].end();
           it != e_it; ++it) {
        ghost_indices.push_back(it->first);
      }
    }

    std::vector< int > offsets(
        static_cast< size_t >(this->nb_subdom_), 0);
    for (size_t i = 1; i < static_cast< size_t >(this->nb_subdom_);
         ++i) {
      offsets[i] = offsets[i - 1] + num_ghost_cells[i - 1];
    }

    std::vector< int > offsets_others(
        static_cast< size_t >(this->nb_subdom_), 0);
    for (size_t i = 1; i < static_cast< size_t >(this->nb_subdom_);
         ++i) {
      offsets_others[i] = offsets_others[i - 1] + num_ghost_cells_others[i - 1];
    }

    // Ghost cell indices which others need from me
    int total_num_ghost_cells_others = 0;
    for (size_t i = 0; i < static_cast< size_t >(this->nb_subdom_);
         ++i) {
      total_num_ghost_cells_others += num_ghost_cells_others[i];
    }

    std::vector< int > ghost_indices_others(
        static_cast< size_t >(total_num_ghost_cells_others), 0);

    MPI_Alltoallv(vec2ptr(ghost_indices), vec2ptr(num_ghost_cells),
                  vec2ptr(offsets), MPI_INT, vec2ptr(ghost_indices_others),
                  vec2ptr(num_ghost_cells_others), vec2ptr(offsets_others),
                  MPI_INT, this->comm_);

    LOG_DEBUG(2, "Exchanged ghost cell indices");

    // Number of ghost dofs for all variables and all ghost cells which
    // I need from others
    std::vector< int > num_dofs_ghost(
        static_cast< size_t >(total_num_ghost_cells * this->nb_fe_), 0);

    // Number of ghost dofs for all variables and all ghost cell which
    // others need from me
    std::vector< int > num_dofs_ghost_others(
        static_cast< size_t >(total_num_ghost_cells_others * this->nb_fe_), 0);

    std::vector< int > offset_ghost_cells_others(
        static_cast< size_t >(this->nb_subdom_), 0);
    for (size_t i = 1; i < static_cast< size_t >(this->nb_subdom_);
         ++i) {
      offset_ghost_cells_others[i] =
          offset_ghost_cells_others[i - 1] + num_ghost_cells_others[i - 1];
    }

    std::vector< int > offset_ghost_cells(
        static_cast< size_t >(this->nb_subdom_), 0);
    for (size_t i = 1; i < static_cast< size_t >(this->nb_subdom_);
         ++i) {
      offset_ghost_cells[i] =
          offset_ghost_cells[i - 1] + num_ghost_cells[i - 1];
    }

    int num_ind_others = 0;
    for (size_t i = 0; i < static_cast< size_t >(this->nb_subdom_);
         ++i) {
      for (size_t k = 0; k < static_cast< size_t >(num_ghost_cells_others[i]);
           ++k) {
        for (size_t l = 0; l < static_cast< size_t >(this->nb_fe_); ++l) {
          num_dofs_ghost_others[(offset_ghost_cells_others[i] + k) * this->nb_fe_ + l] =
              this->nb_dofs_on_cell( l, ghost_indices_others[offset_ghost_cells_others[i] + k]);
          num_ind_others += num_dofs_ghost_others[(offset_ghost_cells_others[i] + k) * this->nb_fe_ + l];
        }
      }
    }

    std::vector< int > num_dofs_ghost_others_vars(
        static_cast< size_t >(this->nb_subdom_), -1);
    for (size_t i = 0; i < static_cast< size_t >(this->nb_subdom_);
         ++i) {
      num_dofs_ghost_others_vars[i] = num_ghost_cells_others[i] * this->nb_fe_;
    }

    for (size_t i = 1; i < static_cast< size_t >(this->nb_subdom_);
         ++i) {
      offsets_others[i] =
          offsets_others[i - 1] + num_dofs_ghost_others_vars[i - 1];
    }

    std::vector< int > num_dofs_ghost_vars(
        static_cast< size_t >(this->nb_subdom_), -1);
    for (size_t i = 0; i < static_cast< size_t >(this->nb_subdom_);
         ++i) {
      num_dofs_ghost_vars[i] = num_ghost_cells[i] * this->nb_fe_;
    }

    for (size_t i = 1; i < static_cast< size_t >(this->nb_subdom_);
         ++i) {
      offsets[i] = offsets[i - 1] + num_dofs_ghost_vars[i - 1];
    }

    MPI_Alltoallv(
        vec2ptr(num_dofs_ghost_others), vec2ptr(num_dofs_ghost_others_vars),
        vec2ptr(offsets_others), MPI_INT, vec2ptr(num_dofs_ghost),
        vec2ptr(num_dofs_ghost_vars), vec2ptr(offsets), MPI_INT, this->comm_);

    LOG_DEBUG(2, "Exchanged number of DoF indices");

    int num_ind = 0;
    for (int i = 0; i < this->nb_subdom_; ++i) {
      for (int k = 0; k < num_ghost_cells[i]; ++k) {
        for (int l = 0; l < this->nb_fe_; ++l) {
          num_ind +=
              num_dofs_ghost[(offset_ghost_cells[i] + k) * this->nb_fe_ + l];
        }
      }
    }

    std::vector< int > recv_num_per_procs(
        static_cast< size_t >(this->nb_subdom_), 0);

    // Prepare final numer_ghost structure
    for (int i = 0; i < this->nb_subdom_; ++i) {
      for (int k = 0; k < num_ghost_cells[i]; ++k) {
        for (int l = 0; l < this->nb_fe_; ++l) {
          numer_ghost[i][ghost_indices[offset_ghost_cells[i] + k]][l].resize(
              num_dofs_ghost[(offset_ghost_cells[i] + k) * this->nb_fe_ + l]);
          recv_num_per_procs[i] +=
              num_dofs_ghost[(offset_ghost_cells[i] + k) * this->nb_fe_ + l];
        }
      }
    }

    LOG_DEBUG(2, "Prepared final numer_ghost structure");

    for (size_t i = 1; i < static_cast< size_t >(this->nb_subdom_);
         ++i) {
      offsets[i] = offsets[i - 1] + recv_num_per_procs[i - 1];
    }

    // Indices which I receive of others
    std::vector< int > recv_indices(static_cast< size_t >(num_ind), 0);

    // Indices which I send to other
    std::vector< int > sent_indices;
    sent_indices.reserve(static_cast< size_t >(num_ind_others));

    std::vector< int > sent_num_per_procs(
        static_cast< size_t >(this->nb_subdom_), 0);

    // Prepare data to send
    for (size_t i = 0; i < this->nb_subdom_; ++i) {
      for (size_t k = 0; k < num_ghost_cells_others[i]; ++k) {
        for (size_t l = 0; l < this->nb_fe_; ++l) {
          std::vector< int > dof_indices;
          this->get_dofs_on_cell(
              l, ghost_indices_others[offset_ghost_cells_others[i] + k],
              dof_indices);
          for (int m = 0; m < dof_indices.size(); ++m) {
            sent_indices.push_back(dof_indices[m]);
          }
          sent_num_per_procs[i] += dof_indices.size();
        }
      }
    }
    assert(sent_indices.size() == num_ind_others);

    LOG_DEBUG(2, "Prepared DoF indices to be sent");

    for (size_t i = 1; i < static_cast< size_t >(this->nb_subdom_);
         ++i) {
      offsets_others[i] = offsets_others[i - 1] + sent_num_per_procs[i - 1];
    }

    MPI_Alltoallv(vec2ptr(sent_indices), vec2ptr(sent_num_per_procs),
                  vec2ptr(offsets_others), MPI_INT, vec2ptr(recv_indices),
                  vec2ptr(recv_num_per_procs), vec2ptr(offsets), MPI_INT,
                  this->comm_);

    LOG_DEBUG(2, "Exchanged DoF indices");

    // Unpack received data
    int ind = 0;
    for (size_t i = 0; i < this->nb_subdom_; ++i) {
      for (size_t k = 0; k < num_ghost_cells[i]; ++k) {
        for (size_t l = 0; l < this->nb_fe_; ++l) {
          for (int m = 0;
               m <
               num_dofs_ghost[(offset_ghost_cells[i] + k) * this->nb_fe_ + l];
               ++m) {
            numer_ghost[i][ghost_indices[offset_ghost_cells[i] + k]][l][m] =
                recv_indices[ind];
            ++ind;
          }
        }
      }
    }
  }

  // First exchange of temporary Dof Ids for ghost layer dofs

  int max_dof_id =
      *(std::max_element(this->numer_cell_2_global_.begin(), this->numer_cell_2_global_.end()));

  std::vector< int > tmp_permutation(static_cast< size_t >(max_dof_id + 1));
  for (size_t i = 0, e_i = tmp_permutation.size(); i != e_i; ++i) {
    tmp_permutation[i] = static_cast< int >(i);
  }

  for (mesh::EntityIterator it = this->mesh_->begin(this->tdim_),
                            e_it = this->mesh_->end(this->tdim_);
       it != e_it; ++it) {
    int subdomain_index;
    it->get("_sub_domain_", &subdomain_index);

    if (subdomain_index != my_subdom_) {
      int ghost_cell_index;
      it->get("_remote_index_", &ghost_cell_index);

      int hostile_tmp_dof_offset = 0;
      for (size_t s = 0, e_s = static_cast< size_t >(subdomain_index); s != e_s;
           ++s) {
        hostile_tmp_dof_offset += ndofs_with_ghost[s];
      }

      for (size_t fe_ind = 0, e_fe_ind = static_cast< size_t >(this->nb_fe_);
           fe_ind != e_fe_ind; ++fe_ind) {
        int size = this->fe_manager_->get_fe(it->index(), fe_ind)->nb_dof_on_cell();

        // Get temporary dof ids from other subdomain
        std::vector< DofID > ghost_layer_dofs(size);
        for (size_t i = 0, e_i = ghost_layer_dofs.size(); i != e_i; ++i) {
          ghost_layer_dofs[i] =
              numer_ghost[subdomain_index][ghost_cell_index][fe_ind][i];
        }

        // Get corresponding dof ids on ghost layer, which need to be updated
        std::vector< DofID > critical_ghost_layer_dofs;
        this->get_dofs_on_cell(fe_ind, it->index(), critical_ghost_layer_dofs);

        for (size_t i = 0, e_i = critical_ghost_layer_dofs.size(); i != e_i;
             ++i) {
          const int cgld_i = critical_ghost_layer_dofs[i];
          if (cgld_i >= tmp_dof_offset + nb_dofs_on_subdom_[my_subdom_] ||
              cgld_i < tmp_dof_offset) {
            const int gld_i = ghost_layer_dofs[i];
            if (gld_i >= hostile_tmp_dof_offset &&
                gld_i <
                    hostile_tmp_dof_offset + nb_dofs_on_subdom_[subdomain_index]) {
              assert(cgld_i >= 0 && cgld_i < tmp_permutation.size());
              tmp_permutation[cgld_i] = gld_i;
            }
          }
        }
      }
    }
  }

  this->apply_permutation(tmp_permutation);

  LOG_DEBUG(2, " Second permutation done ");

  // Update numer field for all subdomains
  numer_ghost.clear();

  // Create subdomain/ghost cell structure of numer_ghost
  for (mesh::EntityIterator it = this->mesh_->begin(this->tdim_),
                            e_it = this->mesh_->end(this->tdim_);
       it != e_it; ++it) {
    int subdomain_index;
    it->get("_sub_domain_", &subdomain_index);

    if (subdomain_index != my_subdom_) {
      // Access/create (on first access) map for remote subdomain
      numer_ghost[subdomain_index];
      int ghost_cell_index;
      it->get("_remote_index_", &ghost_cell_index);
      // Access/create (on first access) map entry for ghost cell
      numer_ghost[subdomain_index][ghost_cell_index];
    }
  }

  {
    // Number of ghost cells which I share with each other process
    std::vector< int > num_ghost_cells(this->nb_subdom_, 0);
    int total_num_ghost_cells = 0;
    for (int i = 0; i < this->nb_subdom_; ++i) {
      num_ghost_cells[i] = numer_ghost[i].size();

      LOG_DEBUG(2, "[" << my_subdom_ << "]: common ghost cells with process "
                       << i << ": " << num_ghost_cells[i]);

      total_num_ghost_cells += numer_ghost[i].size();
    }

    // Number of ghost cells which other processes share with me
    std::vector< int > num_ghost_cells_others(this->nb_subdom_, 0);

    // Exchange number of ghost cells
    MPI_Alltoall(vec2ptr(num_ghost_cells), 1, MPI_INT,
                 vec2ptr(num_ghost_cells_others), 1, MPI_INT, this->comm_);

#ifndef NDEBUG
    for (int i = 0; i < this->nb_subdom_; ++i) {
      LOG_DEBUG(2, "[" << my_subdom_ << "]: process " << i
                       << " requests number of ghost cells: "
                       << num_ghost_cells_others[i]);
    }
#endif

    // Ghost cell indices which I need of others
    std::vector< int > ghost_indices;
    ghost_indices.reserve(total_num_ghost_cells);
    for (int i = 0; i < this->nb_subdom_; ++i) {
      for (typename std::map<
               int, std::map< int, std::vector< int > > >::const_iterator
               it = numer_ghost[i].begin(),
               e_it = numer_ghost[i].end();
           it != e_it; ++it) {
        ghost_indices.push_back(it->first);
      }
    }

    std::vector< int > offsets(this->nb_subdom_, 0);
    for (int i = 1; i < this->nb_subdom_; ++i) {
      offsets[i] = offsets[i - 1] + num_ghost_cells[i - 1];
    }

    std::vector< int > offsets_others(this->nb_subdom_, 0);
    for (int i = 1; i < this->nb_subdom_; ++i) {
      offsets_others[i] = offsets_others[i - 1] + num_ghost_cells_others[i - 1];
    }

    // Ghost cell indices which others need from me
    int total_num_ghost_cells_others = 0;
    for (int i = 0; i < this->nb_subdom_; ++i) {
      total_num_ghost_cells_others += num_ghost_cells_others[i];
    }

    std::vector< int > ghost_indices_others(total_num_ghost_cells_others, 0);

    MPI_Alltoallv(vec2ptr(ghost_indices), vec2ptr(num_ghost_cells),
                  vec2ptr(offsets), MPI_INT, vec2ptr(ghost_indices_others),
                  vec2ptr(num_ghost_cells_others), vec2ptr(offsets_others),
                  MPI_INT, this->comm_);

    LOG_DEBUG(2, "Exchanged ghost cell indices");

    // Number of ghost dofs for all variables and all ghost cells which
    // I need from others
    std::vector< int > num_dofs_ghost(total_num_ghost_cells * this->nb_fe_, 0);

    // Number of ghost dofs for all variables and all ghost cell which
    // others need from me
    std::vector< int > num_dofs_ghost_others(
        total_num_ghost_cells_others * this->nb_fe_, 0);

    std::vector< int > offset_ghost_cells_others(this->nb_subdom_,
                                                 0);
    for (int i = 1; i < this->nb_subdom_; ++i) {
      offset_ghost_cells_others[i] =
          offset_ghost_cells_others[i - 1] + num_ghost_cells_others[i - 1];
    }

    std::vector< int > offset_ghost_cells(this->nb_subdom_, 0);
    for (int i = 1; i < this->nb_subdom_; ++i) {
      offset_ghost_cells[i] =
          offset_ghost_cells[i - 1] + num_ghost_cells[i - 1];
    }

    int num_ind_others = 0;
    for (int i = 0; i < this->nb_subdom_; ++i) {
      for (int k = 0; k < num_ghost_cells_others[i]; ++k) {
        for (int l = 0; l < this->nb_fe_; ++l) {
          num_dofs_ghost_others[(offset_ghost_cells_others[i] + k) *
                                    this->nb_fe_ +
                                l] =
              this->nb_dofs_on_cell(
                  l, ghost_indices_others[offset_ghost_cells_others[i] + k]);
          num_ind_others +=
              num_dofs_ghost_others[(offset_ghost_cells_others[i] + k) *
                                        this->nb_fe_ +
                                    l];
        }
      }
    }

    std::vector< int > num_dofs_ghost_others_vars(this->nb_subdom_,
                                                  -1);
    for (int i = 0; i < this->nb_subdom_; ++i) {
      num_dofs_ghost_others_vars[i] = num_ghost_cells_others[i] * this->nb_fe_;
    }

    for (int i = 1; i < this->nb_subdom_; ++i) {
      offsets_others[i] =
          offsets_others[i - 1] + num_dofs_ghost_others_vars[i - 1];
    }

    std::vector< int > num_dofs_ghost_vars(this->nb_subdom_, -1);
    for (int i = 0; i < this->nb_subdom_; ++i) {
      num_dofs_ghost_vars[i] = num_ghost_cells[i] * this->nb_fe_;
    }

    for (int i = 1; i < this->nb_subdom_; ++i) {
      offsets[i] = offsets[i - 1] + num_dofs_ghost_vars[i - 1];
    }

    MPI_Alltoallv(
        vec2ptr(num_dofs_ghost_others), vec2ptr(num_dofs_ghost_others_vars),
        vec2ptr(offsets_others), MPI_INT, vec2ptr(num_dofs_ghost),
        vec2ptr(num_dofs_ghost_vars), vec2ptr(offsets), MPI_INT, this->comm_);

    LOG_DEBUG(2, "Exchanged number of DoF indices");

    int num_ind = 0;
    for (int i = 0; i < this->nb_subdom_; ++i) {
      for (int k = 0; k < num_ghost_cells[i]; ++k) {
        for (int l = 0; l < this->nb_fe_; ++l) {
          num_ind +=
              num_dofs_ghost[(offset_ghost_cells[i] + k) * this->nb_fe_ + l];
        }
      }
    }

    std::vector< int > recv_num_per_procs(this->nb_subdom_, 0);

    // Prepare final numer_ghost structure
    for (int i = 0; i < this->nb_subdom_; ++i) {
      for (int k = 0; k < num_ghost_cells[i]; ++k) {
        for (int l = 0; l < this->nb_fe_; ++l) {
          numer_ghost[i][ghost_indices[offset_ghost_cells[i] + k]][l].resize(
              num_dofs_ghost[(offset_ghost_cells[i] + k) * this->nb_fe_ + l]);
          recv_num_per_procs[i] +=
              num_dofs_ghost[(offset_ghost_cells[i] + k) * this->nb_fe_ + l];
        }
      }
    }

    LOG_DEBUG(2, "Prepared final numer_ghost structure");

    for (int i = 1; i < this->nb_subdom_; ++i) {
      offsets[i] = offsets[i - 1] + recv_num_per_procs[i - 1];
    }

    // Indices which I receive of others
    std::vector< int > recv_indices(num_ind, 0);

    // Indices which I send to other
    std::vector< int > sent_indices;
    sent_indices.reserve(num_ind_others);

    std::vector< int > sent_num_per_procs(this->nb_subdom_, 0);

    // Prepare data to send
    for (int i = 0; i < this->nb_subdom_; ++i) {
      for (int k = 0; k < num_ghost_cells_others[i]; ++k) {
        for (int l = 0; l < this->nb_fe_; ++l) {
          std::vector< int > dof_indices;
          this->get_dofs_on_cell(
              l, ghost_indices_others[offset_ghost_cells_others[i] + k],
              dof_indices);
          for (int m = 0; m < dof_indices.size(); ++m) {
            sent_indices.push_back(dof_indices[m]);
          }
          sent_num_per_procs[i] += dof_indices.size();
        }
      }
    }
    assert(sent_indices.size() == num_ind_others);

    LOG_DEBUG(2, "Prepared DoF indices to be sent");

    for (int i = 1; i < this->nb_subdom_; ++i) {
      offsets_others[i] = offsets_others[i - 1] + sent_num_per_procs[i - 1];
    }

    MPI_Alltoallv(vec2ptr(sent_indices), vec2ptr(sent_num_per_procs),
                  vec2ptr(offsets_others), MPI_INT, vec2ptr(recv_indices),
                  vec2ptr(recv_num_per_procs), vec2ptr(offsets), MPI_INT,
                  this->comm_);

    LOG_DEBUG(2, "Exchanged DoF indices");

    // Unpack received data
    int ind = 0;
    for (int i = 0; i < this->nb_subdom_; ++i) {
      for (int k = 0; k < num_ghost_cells[i]; ++k) {
        for (int l = 0; l < this->nb_fe_; ++l) {
          for (int m = 0;
               m <
               num_dofs_ghost[(offset_ghost_cells[i] + k) * this->nb_fe_ + l];
               ++m) {
            numer_ghost[i][ghost_indices[offset_ghost_cells[i] + k]][l][m] =
                recv_indices[ind];
            ++ind;
          }
        }
      }
    }
  }

  // Fix temporary Dof Ids on ghost layer to correct dof ids (this step might
  // not always be necessary but is essential to ensure correctness in special
  // cases. See documentation file for more information

  max_dof_id = *(std::max_element(this->numer_cell_2_global_.begin(), this->numer_cell_2_global_.end()));

  std::vector< int > update_permutation(max_dof_id + 1);
  for (size_t i = 0, e_i = update_permutation.size(); i != e_i; ++i) {
    update_permutation[i] = i;
  }

  for (mesh::EntityIterator it = this->mesh_->begin(this->tdim_),
                            e_it = this->mesh_->end(this->tdim_);
       it != e_it; ++it) {
    int subdomain_index;
    it->get("_sub_domain_", &subdomain_index);

    if (subdomain_index != my_subdom_) {
      int ghost_cell_index;
      it->get("_remote_index_", &ghost_cell_index);

      for (size_t fe_ind = 0, e_fe_ind = this->nb_fe_; fe_ind != e_fe_ind; ++fe_ind) {
        // Get dofs from other subdomain
        int size = this->fe_manager_->get_fe(it->index(), fe_ind)->nb_dof_on_cell();

        std::vector< DofID > ghost_layer_dofs(size);
        for (size_t i = 0, e_i = ghost_layer_dofs.size(); i != e_i; ++i) {
          ghost_layer_dofs[i] =
              numer_ghost[subdomain_index][ghost_cell_index][fe_ind][i];
        }

        // Get dofs on ghost layer from view of my subdomain
        std::vector< DofID > critical_ghost_layer_dofs;
        this->get_dofs_on_cell(fe_ind, it->index(), critical_ghost_layer_dofs);

        for (size_t i = 0, e_i = critical_ghost_layer_dofs.size(); i != e_i;
             ++i) {
          const int cgld_i = critical_ghost_layer_dofs[i];
          if (cgld_i >= tmp_dof_offset + nb_dofs_on_subdom_[my_subdom_] ||
              cgld_i < tmp_dof_offset) {
            assert(cgld_i >= 0 && cgld_i < update_permutation.size());
            update_permutation[cgld_i] = ghost_layer_dofs[i];
          }
        }
      }
    }
  }

  this->apply_permutation(update_permutation);

  LOG_DEBUG(2, " Third permutation done ");

  // Finaly calculate real dof_offset and correct numer_ field w.r.t. new offset

  my_dof_offset_ = 0;
  for (size_t s = 0, e_s = my_subdom_; s != e_s; ++s) {
    my_dof_offset_ += nb_dofs_on_subdom_[s];
  }

  std::vector< int > old_dof_offsets(nb_subdom_, 0);
  for (size_t s = 0, e_s = nb_subdom_; s != e_s; ++s) {
    for (size_t t = 0; t != s; ++t) {
      old_dof_offsets[s] += ndofs_with_ghost[t];
    }
  }

  std::vector< int > real_dof_offset(nb_subdom_, 0);
  for (size_t s = 0, e_s = nb_subdom_; s != e_s; ++s) {
    for (size_t t = 0; t != s; ++t) {
      real_dof_offset[s] += nb_dofs_on_subdom_[t];
    }
  }

  int size_of_ownerships = old_dof_offsets[nb_subdom_ - 1] +
                           ndofs_with_ghost[nb_subdom_ - 1];
  std::vector< int > ownerships(size_of_ownerships);

  for (size_t s = 0, e_s = nb_subdom_ - 1; s != e_s; ++s) {
    for (size_t i = old_dof_offsets[s]; i < old_dof_offsets[s + 1]; ++i) {
      ownerships[i] = s;
    }
  }

  for (size_t i = old_dof_offsets[nb_subdom_ - 1],
              e_i = old_dof_offsets[nb_subdom_ - 1] +
                    ndofs_with_ghost[nb_subdom_ - 1];
       i < e_i; ++i) {
    ownerships[i] = nb_subdom_ - 1;
  }

  max_dof_id = *(std::max_element(this->numer_cell_2_global_.begin(), this->numer_cell_2_global_.end()));

  std::vector< int > final_permutation(max_dof_id + 1, -1);
  for (size_t i = 0, e_i = this->numer_cell_2_global_.size(); i != e_i; ++i) {
    int owner = ownerships[this->numer_cell_2_global_[i]];

    if (owner != my_subdom_) {
      final_permutation[this->numer_cell_2_global_[i]] =
          this->numer_cell_2_global_[i] - old_dof_offsets[owner] + real_dof_offset[owner];
    } else {
      final_permutation[this->numer_cell_2_global_[i]] =
          this->numer_cell_2_global_[i] - old_dof_offsets[my_subdom_] + my_dof_offset_;
    }
  }

  this->apply_permutation(final_permutation);

  // Last check if Dof Ids are still greater than -1
  for (size_t i = 0, e_i = this->numer_cell_2_global_.size(); i != e_i; ++i) {
    assert(this->numer_cell_2_global_[i] >= 0);
  }

  // Calculate number of dofs for each variable

  for (size_t fe_ind = 0, e_fe_ind = this->nb_fe_; fe_ind != e_fe_ind; ++fe_ind) {
    int begin_offset = this->numer_cell_2_global_offsets_[fe_ind][0];
    int end_offset;

    if (fe_ind + 1 < this->nb_fe_) {
      end_offset = this->numer_cell_2_global_offsets_[fe_ind + 1][0];
    } else {
      end_offset = this->numer_cell_2_global_.size();
    }

    std::vector< DofID > tmp(end_offset - begin_offset);

    for (size_t i = begin_offset; i < end_offset; ++i) {
      tmp[i - begin_offset] = this->numer_cell_2_global_[i];
    }

    std::sort(tmp.begin(), tmp.end());

    std::vector< DofID >::iterator it = std::unique(tmp.begin(), tmp.end());
    int tmp_size = it - tmp.begin();

    int hostile_dof = 0;
    for (size_t i = 0, e_i = tmp_size; i != e_i; ++i) {
      if (owner_of_dof(tmp[i]) != my_subdom_) {
        hostile_dof++;
      }
    }

    this->local_nb_dofs_for_fe_[fe_ind] -= hostile_dof;
  }

  // consecutive_numbering();
}

template < class DataType, int DIM >
void DofPartitionGlobal< DataType, DIM >::consecutive_numbering() {
  // Fill vector local 2 global and global 2 local map
  local2global_.resize(nb_dofs_incl_ghost_);

  for (size_t fe_ind = 0, e_fe_ind = this->nb_fe_; fe_ind != e_fe_ind; ++fe_ind) {
    int local_dof_cntr = 0;
    // First regular mesh without ghost layer
    for (mesh::EntityIterator it = this->mesh_->begin(this->tdim_),
                              e_it = this->mesh_->end(this->tdim_);
         it != e_it; ++it) {
      int subdomain_index;
      it->get("_sub_domain_", &subdomain_index);

      if (subdomain_index == my_subdom_) {
        std::vector< DofID > global_dofs;
        this->get_dofs_on_cell(fe_ind, it->index(), global_dofs);

        for (size_t i = 0, e_i = global_dofs.size(); i != e_i; ++i) {
          if (global2local_.find(global_dofs[i]) == global2local_.end()) {
            global2local_.insert(
                std::pair< DofID, DofID >(global_dofs[i], local_dof_cntr));
            ++local_dof_cntr;
          }
        }
      }
    }
    // Next: Ghost layer
    for (mesh::EntityIterator it = this->mesh_->begin(this->tdim_),
                              e_it = this->mesh_->end(this->tdim_);
         it != e_it; ++it) {
      int subdomain_index;
      it->get("_sub_domain_", &subdomain_index);

      if (subdomain_index != my_subdom_) {
        std::vector< DofID > global_dofs;
        this->get_dofs_on_cell(fe_ind, it->index(), global_dofs);

        for (size_t i = 0, e_i = global_dofs.size(); i != e_i; ++i) {
          if (global2local_.find(global_dofs[i]) == global2local_.end()) {
            global2local_.insert(
                std::pair< DofID, DofID >(global_dofs[i], local_dof_cntr));
            ++local_dof_cntr;
          }
        }
      }
    }
  }

  // Fill local2global
  std::map< DofID, DofID >::iterator it = global2local_.begin();
  while (it != global2local_.end()) {
    local2global_[(*it).second] = (*it).first;
    ++it;
  }

  // Fill numer_cell_2_local_ field
  this->numer_cell_2_local_.resize(this->numer_cell_2_global_.size(), -1);
  int offset = 0;
  for (size_t fe_ind = 0, e_fe_ind = this->nb_fe_; fe_ind != e_fe_ind; ++fe_ind) 
  {
    for (mesh::EntityIterator it = this->mesh_->begin(this->tdim_),
                              e_it = this->mesh_->end(this->tdim_); it != e_it; ++it) 
    {
      std::vector< DofID > global_dofs;
      this->get_dofs_on_cell(fe_ind, it->index(), global_dofs);

      for (size_t i = 0, e_i = global_dofs.size(); i != e_i; ++i) 
      {
        DofID local_dof;
        global2local(global_dofs[i], &local_dof);
        numer_cell_2_local_[i + offset] = local_dof;
      }

      offset += global_dofs.size();
    }
  }
  
#ifndef NDEBUG
  int diff = 0;
  for (size_t l=0; l<numer_cell_2_local_.size(); ++l)
  {
    diff += std::abs(this->numer_cell_2_local_[l] -  this->numer_cell_2_global_[l]);
  }
  std::cout << " difference numer_sd <> numer_cell_2_subdom : " << diff << std::endl;
  assert (diff == 0);
#endif  
}

template < class DataType, int DIM >
void DofPartitionGlobal< DataType, DIM >::get_dofs_on_cell_local(
    size_t fe_ind, int cell_index, std::vector< DofID > &ids) const {
  // finite element type

  const RefElement< DataType, DIM > &fe_type = *(this->fe_manager_->get_fe(cell_index, fe_ind));

  ids.resize(fe_type.nb_dof_on_cell());

  // loop over DoFs

  for (size_t i = 0, e_i = fe_type.nb_dof_on_cell(); i != e_i; ++i) {
    ids[i] = numer_cell_2_local_[this->numer_cell_2_global_offsets_[fe_ind][cell_index] + i];
  }
}

template < class DataType, int DIM >
void DofPartitionGlobal< DataType, DIM >::apply_permutation_local(
    const std::vector< DofID > &permutation) {

  for (size_t i = 0, e_i = numer_cell_2_local_.size(); i != e_i; ++i) {
    numer_cell_2_local_[i] = permutation[numer_cell_2_local_[i]];
  }

  // Fix local2global and global2local
  std::vector< DofID > l2g_backup(local2global_);
  for (size_t i = 0, e_i = local2global_.size(); i != e_i; ++i) {
    local2global_[i] = l2g_backup[permutation[local2global_[i]]];
  }

  global2local_.clear();
  for (size_t i = 0, e_i = local2global_.size(); i != e_i; ++i) {
    global2local_.insert(std::pair< DofID, DofID >(local2global_[i], i));
  }
}

template class DofPartitionGlobal< double, 3 >;
template class DofPartitionGlobal< double, 2 >;
template class DofPartitionGlobal< double, 1 >;

template class DofPartitionGlobal< float, 3 >;
template class DofPartitionGlobal< float, 2 >;
template class DofPartitionGlobal< float, 1 >;

} // namespace doffem
} // namespace hiflow
