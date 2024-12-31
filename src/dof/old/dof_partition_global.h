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

#ifndef _DOF_PARTITION_GLOBAL_H_
#define _DOF_PARTITION_GLOBAL_H_

#include "common/sorted_array.h"
#include "dof/dof_partition_local.h"
#include <map>
#include <mpi.h>
#include <vector>

namespace hiflow {

namespace mesh {
class Mesh;
}

namespace doffem {

///
/// \class DofPartitionGlobal dof_partition.h
/// \brief This class manages the dof partitioning for parallel computation
/// \author Michael Schick<br>Martin Baumann<br> Simon Gawlok

template < class DataType, int DIM >
class DofPartitionGlobal : public DofPartitionLocal< DataType, DIM > {
public:
  /// Default Constructor
  DofPartitionGlobal();
  /// Destructor
  ~DofPartitionGlobal();

  /// Set MPI Communicator if needed, default value is MPI_COMM_WORLD
  void set_mpi_comm(const MPI_Comm &comm);

  /// Get dof offset of this subdomain
  int my_dof_offset() const;
  
  /// Get subdomain index (subdomain which owns dof)
  int my_subdom() const;
  
  /// Get MPI communicator
  const MPI_Comm &get_mpi_comm() const;

  size_t nb_subdom () const;
  
  /// Returns the number of dofs on a specific subdomain owned by the
  /// subdomain
  size_t nb_dofs_on_subdom(int subdomain) const;
  
  /// Returns the number of dofs on the whole global computational domain
  /// (including the other subdomains)
  size_t nb_dofs_global() const;

  /// Returns dofs on cell for a specific variable. The DofIDs are local
  /// w.r.t. subdomain (Including the dofs which are not owned)
  void get_dofs_on_cell_local(size_t fe_ind, int cell_index, std::vector< DofID > &ids) const;

  /// Create local numbering for dofs on the subdomain
  /// @param[in] order Ordering strategy for DoFs.
  void number_parallel();

  /// Permute the local DofIDs w.r.t. the subdomain by a given permutation
  /// (Including the dofs which are not owned)
  void apply_permutation_local(const std::vector< DofID > &permutation);

  /// Check if global dof is on subdomain (only if dof is owned)
  bool is_dof_on_subdom(DofID global_id) const;
  
  /// Returns the lowest subdomain index, which shares the global dof ID
  int owner_of_dof(DofID global_id) const;

  /// For a given local DofId on a subdomain, this routine computes the global
  /// DofId (including the local dofs which are not owned have local number on
  /// subdomain)
  void local2global(DofID local_id, DofID *global_id) const;
  
  /// For a given global DofId, this routine computes the local DofId on the
  /// current subdomain (including the dofs which are not owned have local
  /// number on subdomain)
  void global2local(DofID global_id, DofID *local_id) const;

private:
  /// Create ownerships of dofs (needed for identification of sceleton dofs)
  void create_ownerships();
  
  /// Renumber dofs according to subdomains to achive global unique numbering
  void renumber();
  
  /// Create local and global correspondences
  void consecutive_numbering();

  /// Subdomain index
  int my_subdom_;
  
  /// Offset for this subdomain
  int my_dof_offset_;
  
  /// Number of dofs including ghost layer
  int nb_dofs_incl_ghost_;
  
  /// Ghost subdomain indices, which are needed in case of point to point
  /// communication
  std::vector< bool > shared_subdomains_;

  /// MPI Communicator
  MPI_Comm comm_;

  /// Total number of subdomains
  int nb_subdom_;

  /// Number of dofs for a specific subdomain owned by the subdomain
  std::vector< int > nb_dofs_on_subdom_;
  
  /// Number of dofs on the global computational domain
  /// (including the other subdomains)
  int gl_nb_dofs_total_;

  /// Vector of ownership
  std::vector< int > ownership_;

  /// Vector storing information for the mapping local 2 global
  std::vector< DofID > local2global_;
  
  /// Map storing information for the mapping global 2 local
  std::map< DofID, DofID > global2local_;
  
  /// Vector storing local (w.r.t. subdomain) numer_ field
  std::vector< DofID > numer_cell_2_local_;
  
};

} // namespace doffem
} // namespace hiflow

#endif
