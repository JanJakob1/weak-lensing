// Copyright (C) 2011-2020 Vincent Heuveline
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

/// \author Aksel Alpay, Martin Wlotzka

#ifndef GMG_RESTRICTION_H
#define GMG_RESTRICTION_H

#include "gmg_interpolation.h"

namespace hiflow {
namespace la {
namespace gmg {

/// Prepares a cell of a vector for restriction by adding contributions from
/// dofs that are ignored during the vector transfer onto the transferred dofs.

template < class LAD, int DIM > 
class LinearCellRestrictor {
public:
  TYPE_FROM_CLASS(LAD, MatrixType);
  TYPE_FROM_CLASS(LAD, VectorType);
  TYPE_FROM_CLASS(LAD, DataType);

  /// Construct object
  /// @param dof_ident A pointer to a dof identification object between the
  /// level on which the vector lives that shall be restricted and the next
  /// coarser level
  /// @param vector A pointer to the vector that shall be prepared for
  /// restriction

  LinearCellRestrictor(
      const boost::shared_ptr< const DofIdentification< LAD, DIM > > &dof_ident,
      VectorType *vector, const std::vector< int > &dirichlet_dofs)
      : dof_ident_(dof_ident), vector_(vector),
        dirichlet_dofs_(dirichlet_dofs) {
    assert(dof_ident);
    assert(vector);
  }

  /// Prepares the restriction matrix.
  /// @param element The cell that shall be used to derive the corresponding
  /// entries in the restriction matrix.
  /// @param matrix The matrix to be prepared.
  void operator()(const Element< DataType, DIM > &element, MatrixType *matrix,
                  std::set< int > &unchanged_dofs_set);

  /// Prepares a cell for restriction
  /// @param element The cell that shall be prepared. It must stem from the
  /// vector that has been supplied in the constructor.
  void operator()(const Element< DataType, DIM > &element);

private:
  /// Checks if a dof on an entity is a dof that must be modified, and if
  /// so, calculates and applies the contributions from the surrounding dofs.
  /// @param processed_dofs A bitmap indicating which dofs on the cell have
  /// already been processed, such that a dof with local id i is located at \c
  /// processed_dofs[i]. If this method modifies the supplied dof, its entry in
  /// the \c processed_dofs bitmap will be set to \c true.
  /// @param dof_partition The dof partition
  /// @param dof_converter A dof converter object for the cell in which the dof
  /// is located
  /// @param dof The global dof id of the dof to be processed
  /// @param var The variable id of the dof
  /// @param ent The (sub)entity on which the dof is located
  inline void process_dof(std::vector< bool > &processed_dofs,
                          const DofPartition< DataType, DIM > &dof_partition,
                          const DofIdConverter< LAD, DIM > &dof_converter,
                          const int dof, const int var, const Entity &ent,
                          MatrixType *matrix);

  /// Checks if a dof on an entity is a dof that must be modified, and if
  /// so, calculates and applies the contributions from the surrounding dofs.
  /// @param processed_dofs A bitmap indicating which dofs on the cell have
  /// already been processed, such that a dof with local id i is located at \c
  /// processed_dofs[i]. If this method modifies the supplied dof, its entry in
  /// the \c processed_dofs bitmap will be set to \c true.
  /// @param dof_partition The dof partition
  /// @param dof_converter A dof converter object for the cell in which the dof
  /// is located
  /// @param dof The global dof id of the dof to be processed
  /// @param var The variable id of the dof
  /// @param ent The (sub)entity on which the dof is located
  inline void process_dof(std::vector< bool > &processed_dofs,
                          const DofPartition< DataType, DIM > &dof_partition,
                          const DofIdConverter< LAD, DIM > &dof_converter,
                          const int dof, const int var, const Entity &ent);

  /// Calculates all contributions to a dof from a given cell.
  /// @param dof_id The global dof id of the dof to be processed
  /// @param var The variable id of the dof
  /// @param cell The cell in which the dof is located
  /// @param dof_converter A dof converter object for this cell
  /// @param weight_sum After a successful call, will contain sum of the weights
  /// of all contributions
  /// @param contributing_dofs The method will append the global ids of the dofs
  /// that contributed to this vector
  /// @param contributing_weights The method will append the weights for
  /// contributing dofs.
  /// @return The sum of the weighted contributions from the surrounding dofs
  inline void get_contributions_from_cell(
      int dof_id, int var, const Element< DataType, DIM > &cell,
      const DofIdConverter< LAD, DIM > &dof_converter,
      std::vector< int > &contributing_dofs,
      std::vector< DataType > &contributing_weights) const;

  /// Calculates all contributions to a dof from a given cell.
  /// @param dof_id The global dof id of the dof to be processed
  /// @param var The variable id of the dof
  /// @param cell The cell in which the dof is located
  /// @param dof_converter A dof converter object for this cell
  /// @param weight_sum After a successful call, will contain sum of the weights
  /// of all contributions
  /// @param contributing_dofs The method will append the global ids of the dofs
  /// that contributed to this vector
  /// @return The sum of the weighted contributions from the surrounding dofs
  inline DataType get_contributions_from_cell(
      int dof_id, int var, const Element< DataType, DIM > &cell,
      const DofIdConverter< LAD, DIM > &dof_converter,
      std::vector< int > &contributing_dofs, DataType &weight_sum) const;

  /// Processes a list of ijk deltas, checks if the corresponding dofs are
  /// located in a given cell, and if so, adds their values and calculates their
  /// weight.
  /// @param delta_list An array of ijk deltas that will be applied onto the
  /// the origin ijk dof id to look for surrounding dofs
  /// @param origin The ijk dof id to which the ijk deltas refer, i.e. the dof
  /// that shall be prepared for restriction
  /// @param var The variable id of the dof
  /// @param dof_converter A dof converter object for this cell
  /// @param weight_of_dof The weight that shall be applied to dofs found via
  /// the delta list
  /// @param contributing_dofs The method will append the global ids of the dofs
  /// that contributed to this vector
  /// @param summed_result_weight After a successful call, will contain the sum
  /// of the weights of the dofs that have contributed
  /// @return The sum of the contributions of the dofs

  template < unsigned N >
  inline void get_contributions_of_delta_list(
      const int (&delta_list)[N][3], const ijk_dof_id &origin, int var,
      const DofIdConverter< LAD, DIM > &dof_converter, const DataType weight_of_dof,
      std::vector< int > &contributing_dofs,
      std::vector< DataType > &contributing_weights) const {
    assert(origin.size() >= 2 && origin.size() <= 3);

    ijk_dof_id current_id = origin;

    for (unsigned i = 0; i < N; ++i) {
      for (unsigned j = 0; j < origin.size(); ++j)
        current_id[j] = origin[j] + delta_list[i][j];

      // check if dof exists and obtain global id
      int local_id = -1;
      if (dof_converter.map_ijk2l(current_id, var, local_id)) {
        int global_id = -1;
        if (dof_converter.map_l2g(local_id, var, global_id)) {
          int lid;
          dof_converter.map_ijk2l(origin, var, lid);
          int gid;
          dof_converter.map_l2g(lid, var, gid);
          //           std::cout << "for dof " << origin[0] << "," << origin[1]
          //           << " with gid " << gid << " found contributing dof " <<
          //           current_id[0] << "," << current_id[1] << " with gid " <<
          //           global_id << " and weight " << weight_of_dof <<
          //           std::endl;
          // if dof has not yet contributed
          if (std::find(contributing_dofs.begin(), contributing_dofs.end(),
                        global_id) == contributing_dofs.end()) {
            //             std::cout << "for dof " << origin[0] << "," <<
            //             origin[1] << " with gid " << gid
            //             << " found contributing dof " << current_id[0] << ","
            //             << current_id[1] << " with gid " << global_id
            //             << " and weight " << weight_of_dof << std::endl;
            contributing_dofs.push_back(global_id);
            contributing_weights.push_back(weight_of_dof);
          }
        }
      }
    }
  }

  /// Processes a list of ijk deltas, checks if the corresponding dofs are
  /// located in a given cell, and if so, adds their values and calculates their
  /// weight.
  /// @param delta_list An array of ijk deltas that will be applied onto the
  /// the origin ijk dof id to look for surrounding dofs
  /// @param origin The ijk dof id to which the ijk deltas refer, i.e. the dof
  /// that shall be prepared for restriction
  /// @param var The variable id of the dof
  /// @param dof_converter A dof converter object for this cell
  /// @param weight_of_dof The weight that shall be applied to dofs found via
  /// the delta list
  /// @param contributing_dofs The method will append the global ids of the dofs
  /// that contributed to this vector
  /// @param summed_result_weight After a successful call, will contain the sum
  /// of the weights of the dofs that have contributed
  /// @return The sum of the contributions of the dofs

  template < unsigned N >
  inline DataType get_contributions_of_delta_list(
      const int (&delta_list)[N][3], const ijk_dof_id &origin, int var,
      const DofIdConverter< LAD, DIM > &dof_converter, const DataType weight_of_dof,
      std::vector< int > &contributing_dofs,
      DataType &summed_result_weight) const {
    // std::cout << "Investigating cell " <<
    // dof_converter.get_element()->get_cell_index() << std::endl;

    assert(origin.size() >= 2 && origin.size() <= 3);

    ijk_dof_id current_id = origin;
    DataType sum_contributions = 0.0;
    summed_result_weight = 0.0;

    for (unsigned i = 0; i < N; ++i) {
      for (unsigned j = 0; j < origin.size(); ++j)
        current_id[j] = origin[j] + delta_list[i][j];

      // check if dof exists and obtain global id
      int local_id = -1;
      if (dof_converter.map_ijk2l(current_id, var, local_id)) {
        int global_id = -1;
        if (dof_converter.map_l2g(local_id, var, global_id)) {
          // if dof has not yet contributed
          if (std::find(contributing_dofs.begin(), contributing_dofs.end(),
                        global_id) == contributing_dofs.end()) {
            // obtain dof value
            DataType value = 0.0;
            vector_->GetValues(&global_id, 1, &value);

            sum_contributions += value * weight_of_dof;
            summed_result_weight += weight_of_dof;

            // std::cout << "Adding " << weight_of_dof << " from " <<
            // dof_converter.get_element()->get_cell_index() << std::endl;

            contributing_dofs.push_back(global_id);
          }
        }
      }
    }

    return sum_contributions;
  }

  /// @return whether a given cell has already been processed
  /// @param current_cell The cell that is currently processed
  /// @param other The cell that shall be checked

  inline bool is_cell_already_processed(const Entity &current_cell,
                                        const Entity &other) const {
    return other.index() < current_cell.index();
  }

  /// @return the number of expected contributions of a dof located in a given
  /// cell
  /// @param cell Specifies the cell type of which in turn is used to determine
  /// expected number of contributions.
  inline unsigned
  get_expected_num_contributions(const Element< DataType, DIM > &cell) const;

  boost::shared_ptr< const DofIdentification< LAD, DIM > > dof_ident_;
  VectorType *vector_;
  OnDemandDofIdConverter< LAD, DIM > on_demand_dof_converter_;
  const std::vector< int > &dirichlet_dofs_;
};

/// Prepares a vector for restriction. The vector will be modified in the
/// process.

template < class LAD, int DIM > 
class LinearRestriction {
public:
  TYPE_FROM_CLASS(LAD, VectorType);

  LinearRestriction< LAD, DIM >() : use_restriction_matrix_(false) {}

  LinearRestriction< LAD, DIM >(const bool use) : use_restriction_matrix_(use) {}

  void use_restriction_matrix(const bool use) { use_restriction_matrix_ = use; }

  template < class ConnectedLevelType >
  void use_transposed_of_interpolation_matrix(ConnectedLevelType &lvl) {
    if (lvl.is_scheduled_to_this_process()) {
      assert(lvl.interpolation_matrix());

      lvl.create_restriction_matrix();
      lvl.restriction_matrix()->CreateTransposedFrom(
          *(lvl.interpolation_matrix().get()));
    }
  }

  template < class ConnectedLevelType >
  void build_restriction_matrix(ConnectedLevelType &lvl) {
    if (lvl.is_scheduled_to_this_process()) {
      lvl.create_restriction_matrix();

      lvl.restriction_matrix()->Zeros();

      ForEachCell< LAD, DIM > for_each_cell(lvl.mesh(), lvl.space());

      LinearCellRestrictor< LAD, DIM > cell_preparator(
          lvl.get_connection_to_next_coarser_grid()->get_dof_identification(),
          lvl.res().get(), lvl.dirichlet_dofs());

      std::set< int > unchanged_dofs_set(lvl.dirichlet_dofs().begin(),
                                         lvl.dirichlet_dofs().end());

      for_each_cell(cell_preparator, lvl.restriction_matrix(),
                    unchanged_dofs_set);

      std::vector< int > unchanged_dofs_vec;

      for (std::set< int >::const_iterator it = unchanged_dofs_set.begin();
           it != unchanged_dofs_set.end(); ++it) {
        if (lvl.space()->dof().is_dof_on_subdom(*it)) {
          unchanged_dofs_vec.push_back(*it);
        }
      }

      lvl.restriction_matrix()->diagonalize_rows(
          vec2ptr(unchanged_dofs_vec), unchanged_dofs_vec.size(), 1.0);

      lvl.restriction_matrix()->Compress();
    }
  }

  /// Run the preparation algorithm.
  /// @param lvl The level in the hierarchy on which the vector is located
  /// @param vector The vector that shall be prepared for restriction. It must
  /// originate from the level given by \c lvl

  template < class ConnectedLevelType >
  void operator()(const ConnectedLevelType &lvl, VectorType &vector,
                  const std::vector< int > &dirichlet_dofs) const {
    if (lvl.is_scheduled_to_this_process()) {
      if (use_restriction_matrix_) {
        assert(lvl.tmp_vector());
        assert(lvl.restriction_matrix());

        lvl.restriction_matrix()->VectorMult(vector, lvl.tmp_vector().get());
        vector.CopyFrom(*(lvl.tmp_vector()));
      } else {
        vector.Update();

        ForEachCell< LAD, DIM > for_each_cell(lvl.mesh(), lvl.space());

        LinearCellRestrictor< LAD, DIM > cell_preparator(
            lvl.get_connection_to_next_coarser_grid()->get_dof_identification(),
            &vector, dirichlet_dofs);

        for_each_cell(cell_preparator);
      }
    }
  }

private:
  bool use_restriction_matrix_;
};

/// Prepares a cell of a vector for restriction by adding contributions from
/// dofs that are ignored during the vector transfer onto the transferred dofs.

template < class LAD, int DIM >
void LinearCellRestrictor< LAD, DIM >::
operator()(const Element< DataType, DIM > &element, MatrixType *matrix,
           std::set< int > &unchanged_dofs_set) {
  unsigned tdim = element.get_cell().tdim();

  const DofIdConverter< LAD, DIM > &dof_converter =
      on_demand_dof_converter_.get_converter(&element);
  const DofPartition< DataType, DIM > &dof_partition = element.get_space().dof();

  std::vector< int > dofs;

  for (int var = 0; var < element.nb_fe(); ++var) {		//TODO: replace var by fe_id?
    std::vector< bool > is_dof_processed(element.nb_dof(var), false);

    // skip dofs on dirichlet boundary
    const int n = dirichlet_dofs_.size();
    for (int i = 0; i < n; ++i) {
      int gid = dirichlet_dofs_[i];
      int lid;
      if (dof_converter.map_g2l(gid, lid)) {
        assert(lid >= 0);
        assert(lid < is_dof_processed.size());
        is_dof_processed[lid] = true;
      }
    }

    // vertices, edges, faces
    for (unsigned dim = 0; dim < tdim; ++dim) {
      int sindex = 0;
      for (IncidentEntityIterator ent = element.get_cell().begin_incident(dim);
           ent != element.get_cell().end_incident(dim); ++ent, ++sindex) {
        // obtain dofs on the subentity
        dof_partition.get_dofs_on_subentity(var, element.cell_index(), dim,
                                            sindex, dofs);

        // Check if the neighboring cells have already processed the dof
        for (IncidentEntityIterator neighbor_cell = ent->begin_incident(tdim);
             neighbor_cell != ent->end_incident(tdim); ++neighbor_cell) {
          if (is_cell_already_processed(element.get_cell(), *neighbor_cell)) {
            // Mark all dofs on this entity as processed
            for (int i = 0; i < dofs.size(); ++i) {
              int current_local_id = -1;
              if (dof_converter.map_g2l(dofs[i], current_local_id)) {
                assert(current_local_id < is_dof_processed.size());
                is_dof_processed[current_local_id] = true;
              }
            }
          }
        }

        for (int i = 0; i < dofs.size(); ++i) {
          int current_local_id = -1;
          if (dof_converter.map_g2l(dofs[i], current_local_id)) {
            if (!is_dof_processed[current_local_id]) {
              process_dof(is_dof_processed, dof_partition, dof_converter,
                          dofs[i], var, *ent, matrix);
            }
          }
        }
      }
    }

    // the cell itself
    element.get_dof_indices(var, dofs);

    for (int i = 0; i < dofs.size(); ++i) {
      // dofs which do not exist on the coarse grid stay unchanged
      if (!dof_ident_->dof_exists_on_coarse_level(dofs[i])) {
        unchanged_dofs_set.insert(dofs[i]);
      }
    }

    for (int i = 0; i < dofs.size(); ++i)
      if (!is_dof_processed[i])
        process_dof(is_dof_processed, dof_partition, dof_converter, dofs[i],
                    var, element.get_cell(), matrix);
  }
}

/// Prepares a cell for restriction
/// @param element The cell that shall be prepared. It must stem from the vector
/// that has been supplied in the constructor.

template < class LAD, int DIM >
void LinearCellRestrictor< LAD, DIM >::
operator()(const Element< DataType, DIM > &element) {
  unsigned tdim = element.get_cell().tdim();

  const DofIdConverter< LAD, DIM > &dof_converter =
      on_demand_dof_converter_.get_converter(&element);
  const DofPartition< DataType, DIM > &dof_partition = element.get_space().dof();

  std::vector< int > dofs;

  for (int var = 0; var < element.nb_fe(); ++var) {		//TODO: replace var by fe_ind
    std::vector< bool > is_dof_processed(element.nb_dof(var), false);

    // skip dofs on dirichlet boundary
    const int n = dirichlet_dofs_.size();
    for (int i = 0; i < n; ++i) {
      int gid = dirichlet_dofs_[i];
      int lid;
      if (dof_converter.map_g2l(gid, lid)) {
        assert(lid >= 0);
        assert(lid < is_dof_processed.size());
        is_dof_processed[lid] = true;
      }
    }

    for (unsigned dim = 0; dim < tdim; ++dim) {
      int sindex = 0;
      for (IncidentEntityIterator ent = element.get_cell().begin_incident(dim);
           ent != element.get_cell().end_incident(dim); ++ent, ++sindex) {
        // obtain dofs
        dof_partition.get_dofs_on_subentity(var, element.cell_index(), dim,
                                            sindex, dofs);

        // Check if the neighboring cells have already processed the dof
        for (IncidentEntityIterator neighbor_cell = ent->begin_incident(tdim);
             neighbor_cell != ent->end_incident(tdim); ++neighbor_cell) {
          if (is_cell_already_processed(element.get_cell(), *neighbor_cell)) {
            // Mark all dofs on this entity as processed
            for (int i = 0; i < dofs.size(); ++i) {
              int current_local_id = -1;
              if (dof_converter.map_g2l(dofs[i], current_local_id)) {
                assert(current_local_id < is_dof_processed.size());
                is_dof_processed[current_local_id] = true;
              }
            }
          }
        }

        for (int i = 0; i < dofs.size(); ++i) {
          int current_local_id = -1;
          if (dof_converter.map_g2l(dofs[i], current_local_id)) {
            if (!is_dof_processed[current_local_id]) {
              process_dof(is_dof_processed, dof_partition, dof_converter,
                          dofs[i], var, *ent);
            }
          }
        }
      }
    }
    element.get_dof_indices(var, dofs);

    for (int i = 0; i < dofs.size(); ++i)
      if (!is_dof_processed[i])
        process_dof(is_dof_processed, dof_partition, dof_converter, dofs[i],
                    var, element.get_cell());
  }
}

/// Checks if a dof on an entity is a dof that must be modified, and if
/// so, calculates and applies the contributions from the surrounding dofs.
/// @param processed_dofs A bitmap indicating which dofs on the cell have
/// already been processed, such that a dof with local id i is located at \c
/// processed_dofs[i]. If this method modifies the supplied dof, its entry in
/// the \c processed_dofs bitmap will be set to \c true.
/// @param dof_partition The dof partition
/// @param dof_converter A dof converter object for the cell in which the dof
/// is located
/// @param dof The global dof id of the dof to be processed
/// @param var The variable id of the dof
/// @param ent The (sub)entity on which the dof is located

template < class LAD, int DIM >
void LinearCellRestrictor< LAD, DIM >::process_dof(
    std::vector< bool > &processed_dofs,
    const DofPartition< DataType, DIM > &dof_partition,
    const DofIdConverter< LAD, DIM > &dof_converter, const int dof, const int var,
    const Entity &ent, MatrixType *matrix) {
  int cell_dim = dof_converter.get_element()->get_cell().tdim();

  int local_id = -1;
  if (dof_converter.map_g2l(dof, local_id)) {
    if (dof_ident_->dof_exists_on_coarse_level(dof)) {
      if (dof_partition.is_dof_on_subdom(dof)) {
        std::vector< int > contributing_dofs;
        std::vector< DataType > contributing_weights;

        contributing_dofs.push_back(dof);
        contributing_weights.push_back(1.0);

        // Get contributions and weights of surrounding dofs
        get_contributions_from_cell(dof, var, *dof_converter.get_element(),
                                    dof_converter, contributing_dofs,
                                    contributing_weights);

        // If we are at a cell border we need to consider the neighboring cells
        if (ent.tdim() < cell_dim) {
          for (IncidentEntityIterator cell = ent.begin_incident(cell_dim);
               cell != ent.end_incident(cell_dim); ++cell) {
            if (cell->index() !=
                dof_converter.get_element()->cell_index()) {
              Element< DataType, DIM > neighbor_cell(
                  dof_converter.get_element()->get_space(), cell->index());

              DofIdConverter< LAD, DIM > neighbor_converter(&neighbor_cell);

              get_contributions_from_cell(dof, var, neighbor_cell,
                                          neighbor_converter, contributing_dofs,
                                          contributing_weights);
            }
          }
        }

        std::map< int, DataType > sorted;
        for (int i = 0; i < contributing_dofs.size(); ++i)
          sorted.insert(
              std::make_pair(contributing_dofs[i], contributing_weights[i]));

        contributing_dofs.clear();
        contributing_weights.clear();

        for (typename std::map< int, DataType >::const_iterator it =
                 sorted.begin();
             it != sorted.end(); ++it) {
          contributing_dofs.push_back(it->first);
          contributing_weights.push_back(it->second);
        }
        matrix->Add(&dof, 1, vec2ptr(contributing_dofs),
                    contributing_dofs.size(), vec2ptr(contributing_weights));
      }
    }
    processed_dofs[local_id] = true;
  }
}

/// Checks if a dof on an entity is a dof that must be modified, and if
/// so, calculates and applies the contributions from the surrounding dofs.
/// @param processed_dofs A bitmap indicating which dofs on the cell have
/// already been processed, such that a dof with local id i is located at \c
/// processed_dofs[i]. If this method modifies the supplied dof, its entry in
/// the \c processed_dofs bitmap will be set to \c true.
/// @param dof_partition The dof partition
/// @param dof_converter A dof converter object for the cell in which the dof
/// is located
/// @param dof The global dof id of the dof to be processed
/// @param var The variable id of the dof
/// @param ent The (sub)entity on which the dof is located

template < class LAD, int DIM >
void LinearCellRestrictor< LAD, DIM >::process_dof(
    std::vector< bool > &processed_dofs,
    const DofPartition< DataType, DIM > &dof_partition,
    const DofIdConverter< LAD, DIM > &dof_converter, const int dof, const int var,
    const Entity &ent) {
  int cell_dim = dof_converter.get_element()->get_cell().tdim();

  int local_id = -1;
  if (dof_converter.map_g2l(dof, local_id)) {
    if (dof_ident_->dof_exists_on_coarse_level(dof)) {
      if (dof_partition.is_dof_on_subdom(dof)) {
        std::vector< int > contributing_dofs;

        DataType contributions = 0.0;
        DataType weight_sum = 0.0;

        // Get contributions and weights of surrounding dofs
        contributions += get_contributions_from_cell(
            dof, var, *dof_converter.get_element(), dof_converter,
            contributing_dofs, weight_sum);

        // If we are at a cell border we need to consider the neighboring cells
        if (ent.tdim() < cell_dim) {
          for (IncidentEntityIterator cell = ent.begin_incident(cell_dim);
               cell != ent.end_incident(cell_dim); ++cell) {
            if (cell->index() !=
                dof_converter.get_element()->cell_index()) {
              Element< DataType, DIM > neighbor_cell(
                  dof_converter.get_element()->get_space(), cell->index());

              DofIdConverter< LAD, DIM > neighbor_converter(&neighbor_cell);

              DataType weights = 0.0;
              contributions += get_contributions_from_cell(
                  dof, var, neighbor_cell, neighbor_converter,
                  contributing_dofs, weights);
              weight_sum += weights;
            }
          }
        }

        // Add previous value
        DataType previous_value = 0.0;
        vector_->GetValues(&dof, 1, &previous_value);

        contributions += previous_value;
        weight_sum += 1.0;

        // std::cout << "dof = " << dof << " prev = " << previous_value << "
        // contributions = " << contributions << " weight = " << weight_sum <<
        // std::endl;

        DataType new_value = contributions;

        vector_->SetValues(&dof, 1, &new_value);
        // std::cout << "---------------------\n";
      }
    }
    processed_dofs[local_id] = true;
  }
}

/// Calculates all contributions to a dof from a given cell.
/// @param dof_id The global dof id of the dof to be processed
/// @param var The variable id of the dof
/// @param cell The cell in which the dof is located
/// @param dof_converter A dof converter object for this cell
/// @param weight_sum After a successful call, will contain sum of the weights
/// of all contributions
/// @param contributing_dofs The method will append the global ids of the dofs
/// that contributed to this vector
/// @param contributing_weights The method will append the weights for
/// contributing dofs.
/// @return The sum of the weighted contributions from the surrounding dofs

template < class LAD, int DIM >
void LinearCellRestrictor< LAD, DIM >::get_contributions_from_cell(
    int dof_id, int var, const Element< DataType, DIM > &cell,
    const DofIdConverter< LAD, DIM > &dof_converter,
    std::vector< int > &contributing_dofs,
    std::vector< DataType > &contributing_weights) const {
  int local_id = -1;
  if (dof_converter.map_g2l(dof_id, local_id)) {
    ijk_dof_id ijk_id;
    if (dof_converter.map_l2ijk(local_id, var, ijk_id)) {
      CellType::Tag t = cell.get_cell().cell_type().tag();

      switch (t) {
      case CellType::QUADRILATERAL:
        get_contributions_of_delta_list(
            NearestKnownDofsFinder< LAD, DIM >::QDeltas::deltas2d_on_axis, ijk_id,
            var, dof_converter, 0.5, contributing_dofs, contributing_weights);

        get_contributions_of_delta_list(
            NearestKnownDofsFinder< LAD, DIM >::QDeltas::deltas2d_diagonals, ijk_id,
            var, dof_converter, 0.25, contributing_dofs, contributing_weights);
        break;
      case CellType::TRIANGLE:
        get_contributions_of_delta_list(
            NearestKnownDofsFinder< LAD, DIM >::PDeltas::deltas2d, ijk_id, var,
            dof_converter, 0.5, contributing_dofs, contributing_weights);

        break;
      case CellType::HEXAHEDRON:
        get_contributions_of_delta_list(
            NearestKnownDofsFinder< LAD, DIM >::QDeltas::deltas3d_on_axis, ijk_id,
            var, dof_converter, 0.5, contributing_dofs, contributing_weights);

        get_contributions_of_delta_list(
            NearestKnownDofsFinder< LAD, DIM >::QDeltas::deltas3d_two_axis_diagonals,
            ijk_id, var, dof_converter, 0.25, contributing_dofs,
            contributing_weights);

        get_contributions_of_delta_list(
            NearestKnownDofsFinder<
                LAD, DIM >::QDeltas::deltas3d_three_axis_diagonals,
            ijk_id, var, dof_converter, 0.125, contributing_dofs,
            contributing_weights);
        break;
      case CellType::TETRAHEDRON:
        // TODO
        throw std::runtime_error(
            "Restriction Preparator: Tetrahedrons not yet supported.");
        break;
      case CellType::PYRAMID:
        // TODO
        throw std::runtime_error(
            "Restriction Preparator: Pyramids not yet supported.");
        break;
      case CellType::POINT:
        // TODO
        throw std::runtime_error(
            "Restriction Preparator: Points not yet supported.");
        break;
      case CellType::LINE:
        // TODO
        throw std::runtime_error(
            "Restriction Preparator: Lines not yet supported.");
        break;
      default:
        throw std::runtime_error("Not a implemented element.");
        break;
      }
    }
  }
}

/// Calculates all contributions to a dof from a given cell.
/// @param dof_id The global dof id of the dof to be processed
/// @param var The variable id of the dof
/// @param cell The cell in which the dof is located
/// @param dof_converter A dof converter object for this cell
/// @param weight_sum After a successful call, will contain sum of the weights
/// of all contributions
/// @param contributing_dofs The method will append the global ids of the dofs
/// that contributed to this vector
/// @return The sum of the weighted contributions from the surrounding dofs

template < class LAD, int DIM >
typename LinearCellRestrictor< LAD, DIM >::DataType
LinearCellRestrictor< LAD, DIM >::get_contributions_from_cell(
    int dof_id, int var, const Element< DataType, DIM > &cell,
    const DofIdConverter< LAD, DIM > &dof_converter,
    std::vector< int > &contributing_dofs, DataType &weight_sum) const {
  weight_sum = 0.0;

  DataType result = 0.0;

  int local_id = -1;
  if (dof_converter.map_g2l(dof_id, local_id)) {
    ijk_dof_id ijk_id;
    if (dof_converter.map_l2ijk(local_id, var, ijk_id)) {
      CellType::Tag t = cell.get_cell().cell_type().tag();

      DataType weight_temp = 0.0;
      switch (t) {
      case CellType::QUADRILATERAL:
        result += get_contributions_of_delta_list(
            NearestKnownDofsFinder< LAD, DIM >::QDeltas::deltas2d_on_axis, ijk_id,
            var, dof_converter, 0.5, contributing_dofs, weight_temp);
        weight_sum += weight_temp;

        result += get_contributions_of_delta_list(
            NearestKnownDofsFinder< LAD, DIM >::QDeltas::deltas2d_diagonals, ijk_id,
            var, dof_converter, 0.25, contributing_dofs, weight_temp);
        weight_sum += weight_temp;

        break;
      case CellType::TRIANGLE:
        result += get_contributions_of_delta_list(
            NearestKnownDofsFinder< LAD, DIM >::PDeltas::deltas2d, ijk_id, var,
            dof_converter, 0.5, contributing_dofs, weight_sum);

        break;
      case CellType::HEXAHEDRON:
        result += get_contributions_of_delta_list(
            NearestKnownDofsFinder< LAD, DIM >::QDeltas::deltas3d_on_axis, ijk_id,
            var, dof_converter, 0.5, contributing_dofs, weight_temp);
        weight_sum += weight_temp;

        result += get_contributions_of_delta_list(
            NearestKnownDofsFinder< LAD, DIM >::QDeltas::deltas3d_two_axis_diagonals,
            ijk_id, var, dof_converter, 0.25, contributing_dofs, weight_temp);
        weight_sum += weight_temp;

        result += get_contributions_of_delta_list(
            NearestKnownDofsFinder<
                LAD, DIM >::QDeltas::deltas3d_three_axis_diagonals,
            ijk_id, var, dof_converter, 0.125, contributing_dofs, weight_temp);
        weight_sum += weight_temp;

        break;
      case CellType::TETRAHEDRON:
        // TODO
        throw std::runtime_error(
            "Restriction Preparator: Tetrahedrons not yet supported.");
        break;
      case CellType::PYRAMID:
        // TODO
        throw std::runtime_error(
            "Restriction Preparator: Pyramids not yet supported.");
        break;
      case CellType::POINT:
        // TODO
        throw std::runtime_error(
            "Restriction Preparator: Points not yet supported.");
        break;
      case CellType::LINE:
        // TODO
        throw std::runtime_error(
            "Restriction Preparator: Lines not yet supported.");
        break;
      default:
        throw std::runtime_error("Not a implemented element.");
        break;
      }

      return result;
    }
  }

  return 0.0;
}

/// @return the number of expected contributions of a dof located in a given
/// cell
/// @param cell Specifies the cell type of which in turn is used to determine
/// expected number of contributions.

template < class LAD, int DIM >
unsigned LinearCellRestrictor< LAD, DIM >::get_expected_num_contributions(
    const Element< DataType, DIM > &cell) const {
  CellType::Tag t = cell.get_cell().cell_type().tag();

  switch (t) {
  case CellType::QUADRILATERAL:
    return 8;
    break;
  case CellType::TRIANGLE:
    return 6;
    break;
  case CellType::HEXAHEDRON:
    return 26;
    break;
  case CellType::TETRAHEDRON:
    return 12;
    break;
  default:
    throw std::runtime_error("Elements are not supported yet.");
    break;
  }
}

} // namespace gmg
} // namespace la
} // namespace hiflow

#endif
