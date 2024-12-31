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

#ifndef _DOF_NUMBERING_STRATEGY_H_
#define _DOF_NUMBERING_STRATEGY_H_

#include <map>
#include <vector>

#include "common/vector_algebra.h"
#include "dof/dof_fem_types.h"
#include "dof/dof_interpolation.h"
#include "dof/dof_interpolation_pattern.h"
#include "dof/fe_interface_pattern.h"
#include "dof/degree_of_freedom.h"
#include "fem/femanager.h"
#include "mesh/entity.h"
#include "mesh/mesh.h"

namespace hiflow {
namespace doffem {

template < class DataType, int DIM > class FEManager;

template < class DataType, int DIM > class FEType;

template < class DataType, int DIM > class DegreeOfFreedom;

/// Enumeration of different DoF ordering strategies. HIFLOW_CLASSIC refers to
/// the DoF numbering as always done in HiFlow3. The other two options allow
/// permutations of the classic numbering by means of the Cuthill-McKee and the
/// King method, respectively.

enum DOF_ORDERING { HIFLOW_CLASSIC, CUTHILL_MCKEE, KING };

/// \author Michael Schick<br>Martin Baumann<br>Simon Gawlok

template < class DataType, int DIM > 
class NumberingStrategy {
public:
  typedef Vec<DIM, DataType> Coord;

  /// Constructor

  NumberingStrategy() {}
  /// Destructor

  virtual ~NumberingStrategy() {}

  /// Initialization. Here, the critical member variables of DegreeOfFreedom are
  /// being given to NumberingStrategy, such that the implementation class of
  /// the function void number() can calculate all neccessary information and
  /// store it in these variables
  void initialize(DegreeOfFreedom< DataType, DIM > &dof);

  /// Kernel of numbering strategy. Here the user can specify in an inherited
  /// class his wishes for some numbering procedure, when dealing for example
  /// with varying finite element types
  /// @param[in] order Ordering strategy for DoFs.
  virtual void number(DOF_ORDERING order = HIFLOW_CLASSIC) = 0;

  /// Helper function which permutes data within the interpolation container and
  /// numer_ field
  void apply_permutation(const std::vector< DofID > &permutation,
                         const std::string & = "");

protected:
  /// Topological dimension
  int tdim_;
  /// Total number of variables
  int dof_nvar_;

  /// The DegreeOfFreedom class used throughout the numeration procedure
  DegreeOfFreedom< DataType, DIM > *dof_;

  /// Const pointer to mesh
  const mesh::Mesh *mesh_;

  /// FEManager on the given mesh
  FEManager< DataType, DIM > const *fe_manager_;

  /// Holds DoF IDs, needed for mapl2g
  std::vector< DofID > *numer_;

  /// Offset container for numer_, needed for mapl2g
  std::vector< std::vector< int > > *numer_offsets_cell_dofvarloc_;

  /// Interpolation Container, which stores the interpolation weights
  DofInterpolation *dof_interpolation_;

  /// Update number_of_dofs_total_ and number_of_dofs_for_var_
  void update_number_of_dofs(const std::string & = "");

  /// Total number of dofs for all variables
  int number_of_dofs_total_;

  /// Total number of dofs per variable
  std::vector< int > number_of_dofs_for_var_;

  /// Get vertex points of mesh entity.
  void get_points(const mesh::Entity &entity, std::vector< Coord > &points);
};

template < class DataType, int DIM >
void NumberingStrategy< DataType, DIM >::initialize(DegreeOfFreedom< DataType, DIM > &dof) {
  dof_ = &dof;
  mesh_ = &(dof_->get_mesh());
  fe_manager_ = &(dof_->get_fe_manager());

  numer_ = dof.numer();
  numer_offsets_cell_dofvarloc_ = dof.numer_offsets_cell_dofvarloc();
  dof_interpolation_ = &(dof.dof_interpolation());

  tdim_ = mesh_->tdim();
  dof_nvar_ = fe_manager_->aaa_get_nb_var();
}

/// \param description is an optional parameter that should describe what
///                    the permutation represents

template < class DataType, int DIM >
void NumberingStrategy< DataType, DIM >::apply_permutation(
    const std::vector< DofID > &permutation, const std::string &description) {
  // apply permutation to mapl2g

  // DoF IDs are used in numer_ only

  for (size_t i = 0, e_i = numer_->size(); i != e_i; ++i) {
    (*numer_)[i] = permutation[(*numer_)[i]];
  }

  // apply permutation to DofInterpolation

  dof_interpolation_->apply_permutation(permutation);

  // calculate number of dofs, as this could have changed

  update_number_of_dofs(description);
}

/// \param description is an optional parameter that should describe what
///                    the permutation represents

template < class DataType, int DIM >
void NumberingStrategy< DataType, DIM >::update_number_of_dofs(
    const std::string &description) {
  // Calculate number of DoFs

  number_of_dofs_total_ = 0;
  for (size_t i = 0, e_i = numer_->size(); i != e_i; ++i) {
    if ((*numer_)[i] > number_of_dofs_total_) {
      number_of_dofs_total_ = (*numer_)[i];
    }
  }
  ++number_of_dofs_total_;

  // Calculate number of Dofs for each variable

  number_of_dofs_for_var_.resize(dof_nvar_, 0);
  for (size_t var = 0; var != dof_nvar_; ++var) {
    int begin_offset = (*numer_offsets_cell_dofvarloc_)[var][0];
    int end_offset;

    if (var < dof_nvar_ - 1) {
      end_offset = (*numer_offsets_cell_dofvarloc_)[var + 1][0];
    } else {
      end_offset = numer_->size();
    }

    for (size_t i = begin_offset; i < end_offset; ++i) {
      if ((*numer_)[i] > number_of_dofs_for_var_[var]) {
        number_of_dofs_for_var_[var] = (*numer_)[i];
      }
    }

    for (size_t j = 0; j != var; ++j) {
      number_of_dofs_for_var_[var] -= number_of_dofs_for_var_[j];
    }

    ++number_of_dofs_for_var_[var];
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
void NumberingStrategy< DataType, DIM >::get_points(const mesh::Entity &entity,
                                                    std::vector< Coord > &points) {
  points.reserve(points.size() + entity.num_vertices());
  for (size_t p = 0; p != entity.num_vertices(); ++p) {
    std::vector< DataType > coords;
    entity.get_coordinates(coords);
    points.push_back( Vec<DIM, DataType>(coords) );
  }
}

} // namespace doffem
} // namespace hiflow

#endif
