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

/// \author Staffan Ronnas, Simon Gawlok

#ifndef _STANDARD_ASSEMBLY_H_
#define _STANDARD_ASSEMBLY_H_

#include <vector>
#include "assembly/assembly_utils.h"
#include "assembly/global_assembler.h"
#include "assembly/generic_assembly_algorithm.h"
#include "common/pointers.h"
#include "common/sort_permutation.h"
#include "mesh/attributes.h"
#include "mesh/types.h"
#include "mesh/iterator.h"
#include "space/element.h"
#include "space/vector_space.h"

namespace hiflow {

//////////////// StandardAssembly helper functions ////////////////

template < template < class, class, int > class AlgorithmType, class DataType, int DIM >
class StandardScalarAssembly
    : public AssemblyAlgorithmBase<AlgorithmType, StandardScalarAssembly< AlgorithmType, DataType, DIM >, DataType, DIM > 
{
public:
  typedef DataType LocalObjectType;
  typedef hiflow::Quadrature< DataType > QuadratureType;
  typedef AssemblyAlgorithmBase< AlgorithmType, StandardScalarAssembly< AlgorithmType, DataType, DIM >, DataType, DIM > Base;

  // Name resolution does not manage to get the base members, so
  // we must do it ourselves.
  using Base::curr_;
  using Base::elem_;
  using Base::has_next;
  using Base::space_;
  using Base::traversal_;

  StandardScalarAssembly(const VectorSpace< DataType, DIM > &space,
                         std::vector< DataType > &values)
      : Base(space), values_(values) {
    const size_t num_elements = this->traversal_.size();
    this->values_.resize(num_elements, 0.);

    this->remove_non_local_elements();
    sort_elements(space, this->traversal_);
  }

  /// The has_next() and next() functions are overloaded in
  /// order to skip elements that do not belong to the local
  /// subdomain. This is done by setting the corresponding
  /// entries in the traversal_ array to -1, and later skipping
  /// those items.

  const Element< DataType, DIM > &next() {
    assert(this->has_next());

    this->elem_ = Element< DataType, DIM >(this->space_, this->traversal_[this->curr_]);

    ++(this->curr_);

    return this->elem_;
  }

  void add(const Element< DataType, DIM > &element,
           const LocalObjectType &local_val) {
    this->values_[element.cell_index()] += local_val;
  }

private:
  /// Remove non_local elements

  void remove_non_local_elements() {

    const mesh::Mesh &mesh = this->space_.mesh();

    if (!mesh.has_attribute("_remote_index_", mesh.tdim())) {
      // If the "_remote_index_" attribute does not exist, we
      // assume that there are no ghost cells.
      return;
    }

    int remote_index;
    int index;
    std::vector< int > traversal_tmp;
    traversal_tmp.reserve(mesh.num_entities(mesh.tdim()));

    for (mesh::EntityIterator it_cell = mesh.begin(mesh.tdim()),
                              e_it_cell = mesh.end(mesh.tdim());
         it_cell != e_it_cell; ++it_cell) {
      // test if cell on subdomain
      it_cell->get("_remote_index_", &remote_index);
      if (remote_index == -1) {
        index = it_cell->index();
        traversal_tmp.push_back(index);
      }
    }
    this->traversal_ = traversal_tmp;
  }

  std::vector< DataType > &values_;
};

template < template < class, class, int > class AlgorithmType, class DataType, int DIM >
class StandardMultipleScalarAssembly
    : public AssemblyAlgorithmBase<AlgorithmType,
                                   StandardMultipleScalarAssembly< AlgorithmType, DataType, DIM >,
                                   DataType, DIM > 
{
public:
  typedef std::vector< DataType > LocalObjectType;
  typedef hiflow::Quadrature< DataType > QuadratureType;
  typedef AssemblyAlgorithmBase< AlgorithmType, StandardMultipleScalarAssembly< AlgorithmType, DataType, DIM >, DataType, DIM > Base;

  // Name resolution does not manage to get the base members, so
  // we must do it ourselves.
  using Base::curr_;
  using Base::elem_;
  using Base::has_next;
  using Base::space_;
  using Base::traversal_;

  StandardMultipleScalarAssembly(const VectorSpace< DataType, DIM > &space,
                                 std::vector< std::vector< DataType > > &values,
                                 const size_t num_scalars)
      : Base(space), num_scalars_(num_scalars), values_(values) {
    const size_t num_elements = this->traversal_.size();
    this->values_.resize(num_elements);
    for (size_t l = 0; l < num_elements; ++l) {
      this->values_[l].resize(this->num_scalars_, 0.);
    }
    this->remove_non_local_elements();
    sort_elements(space, this->traversal_);
  }

  /// The has_next() and next() functions are overloaded in
  /// order to skip elements that do not belong to the local
  /// subdomain. This is done by setting the corresponding
  /// entries in the traversal_ array to -1, and later skipping
  /// those items.

  const Element< DataType, DIM > &next() {
    assert(this->has_next());

    this->elem_ = Element< DataType, DIM >(this->space_, this->traversal_[this->curr_]);

    ++(this->curr_);

    return this->elem_;
  }

  void add(const Element< DataType, DIM > &element, const LocalObjectType &local_val) 
  {
    for (size_t l = 0; l < this->num_scalars_; ++l) 
    {
      this->values_[element.cell_index()][l] += local_val[l];
    }
  }

private:
  /// Remove non_local elements

  void remove_non_local_elements() {

    const mesh::Mesh &mesh = this->space_.mesh();

    if (!mesh.has_attribute("_remote_index_", mesh.tdim())) {
      // If the "_remote_index_" attribute does not exist, we
      // assume that there are no ghost cells.
      return;
    }

    int remote_index;
    int index;
    std::vector< int > traversal_tmp;
    traversal_tmp.reserve(mesh.num_entities(mesh.tdim()));

    for (mesh::EntityIterator it_cell = mesh.begin(mesh.tdim()),
                              e_it_cell = mesh.end(mesh.tdim());
         it_cell != e_it_cell; ++it_cell) {
      // test if cell on subdomain
      it_cell->get("_remote_index_", &remote_index);
      if (remote_index == -1) {
        index = it_cell->index();
        traversal_tmp.push_back(index);
      }
    }
    this->traversal_ = traversal_tmp;
  }
  size_t num_scalars_;
  std::vector< std::vector< DataType > > &values_;
};

template < template < class, class, int > class AlgorithmType, class DataType, int DIM >
class StandardVectorAssembly
    : public AssemblyAlgorithmBase<AlgorithmType, StandardVectorAssembly< AlgorithmType, DataType, DIM >, DataType, DIM > 
{
public:
  typedef std::vector< DataType > LocalObjectType;
  typedef hiflow::Quadrature< DataType > QuadratureType;
  typedef AssemblyAlgorithmBase< AlgorithmType, StandardVectorAssembly< AlgorithmType, DataType, DIM >, DataType, DIM > Base;

  // Name resolution does not manage to get the base members, so
  // we must do it ourselves.
  using Base::dof_;
  using Base::space_;
  using Base::traversal_;

  StandardVectorAssembly(
      const VectorSpace< DataType, DIM > &space,
      typename GlobalAssembler< DataType, DIM >::GlobalVector &vec)
      : Base(space), vector_(vec) {

    sort_elements(space, this->traversal_);
  }

  void add(const Element< DataType, DIM > &element,
           const LocalObjectType &local_vec) {

    assert (!contains_nan(local_vec));
    
    const size_t num_dofs = this->dof_.size();

    std::vector< int > dofs_sort_permutation;

    // get permutation for sorting dofs
    sortingPermutation(this->dof_, dofs_sort_permutation);

    // create row array
    std::vector< int > row_indices;
    row_indices.reserve(num_dofs);

    LocalObjectType local_vec_sorted;
    local_vec_sorted.reserve(num_dofs);

    std::vector< DataType > dof_factors;
    this->space_.dof().get_dof_factors_on_cell (element.cell_index(), dof_factors);
    assert (dof_factors.size() == num_dofs);
    
    for (size_t i = 0; i != num_dofs; ++i) 
    {
      const int dof_sort_perm = dofs_sort_permutation[i];
      const int dof_ind = this->dof_[dof_sort_perm];
      if (this->space_.dof().is_dof_on_subdom(dof_ind)) 
      {
        row_indices.push_back(dof_ind);
        local_vec_sorted.push_back(dof_factors[dof_sort_perm] * local_vec[dof_sort_perm]);
      }
    }
            
    // Add local to global vector
    if (!row_indices.empty()) {
      this->vector_.Add(vec2ptr(row_indices), row_indices.size(),
                        vec2ptr(local_vec_sorted));
    }
  }

private:
  typename GlobalAssembler< DataType, DIM >::GlobalVector &vector_;
};

template < template < class, class, int > class AlgorithmType, class DataType, int DIM >
class StandardMatrixAssembly
    : public AssemblyAlgorithmBase<AlgorithmType, StandardMatrixAssembly< AlgorithmType, DataType, DIM >, DataType, DIM > 
{
public:
  typedef la::SeqDenseMatrix< DataType > LocalObjectType;
  typedef Quadrature< DataType > QuadratureType;
  typedef AssemblyAlgorithmBase<AlgorithmType, StandardMatrixAssembly< AlgorithmType, DataType, DIM >, DataType, DIM > Base;

  // Name resolution does not manage to get the base members, so
  // we must do it ourselves.
  using Base::dof_;
  using Base::space_;
  using Base::traversal_;

  StandardMatrixAssembly(
      const VectorSpace< DataType, DIM > &space,
      typename GlobalAssembler< DataType, DIM >::GlobalMatrix &matrix)
      : Base(space), matrix_(matrix) {
    sort_elements(space, this->traversal_);
  }

  void add(const Element< DataType, DIM > &element,
           const LocalObjectType &local_mat) {
    const size_t num_dofs = this->dof_.size();

    std::vector< int > dofs_sort_permutation;
    std::vector< int > dofs_sorted(num_dofs);

    // get permutation for sorting dofs
    sortingPermutation(this->dof_, dofs_sort_permutation);

    // fill sorted dof array
    for (size_t i = 0; i != num_dofs; ++i) {
      dofs_sorted[i] = this->dof_[dofs_sort_permutation[i]];
    }

    std::vector< DataType > dof_factors;
    this->space_.dof().get_dof_factors_on_cell (element.cell_index(), dof_factors);
    assert (dof_factors.size() == num_dofs);
    
    // create row array
    std::vector< int > row_indices;
    std::vector< int > row_permutation;
    row_indices.reserve(num_dofs);
    row_permutation.reserve(num_dofs);
    for (size_t i = 0; i != num_dofs; ++i) 
    {
      const int dof_sort_perm = dofs_sort_permutation[i];
      const int dof_ind = this->dof_[dof_sort_perm];
      if (this->space_.dof().is_dof_on_subdom(dof_ind)) 
      {
        row_indices.push_back(dof_ind);
        row_permutation.push_back(dof_sort_perm);
      }
    }

    // fill reduced and sorted local matrix
    // TODO: make local_mat_sorted_reduced as class member for performance issues
    LocalObjectType local_mat_sorted_reduced;
    if (!row_indices.empty() && num_dofs > 0) 
    {
      local_mat_sorted_reduced.Resize(row_indices.size(), num_dofs);
      for (size_t i = 0, i_e = row_indices.size(); i != i_e; ++i) 
      {
        const int row_ind = row_permutation[i];
//std::cout << " CG " << row_ind << " <-> " << i << " <-> " << dofs_sort_permutation[i] << std::endl;
        for (size_t j = 0, j_e = num_dofs; j != j_e; ++j) 
        {
          const int col_ind = dofs_sort_permutation[j];
          local_mat_sorted_reduced(i, j) =  dof_factors[row_ind] 
                                          * dof_factors[col_ind] 
                                          * local_mat(row_ind, col_ind);
        }
      }

      // Add local to global matrix
      assert (row_indices.size() > 0);
      assert (dofs_sorted.size() > 0);
      
      this->matrix_.Add(vec2ptr(row_indices), row_indices.size(),
                        vec2ptr(dofs_sorted), dofs_sorted.size(),
                        &local_mat_sorted_reduced(0, 0));
    }
  }

private:
  typename GlobalAssembler< DataType, DIM >::GlobalMatrix &matrix_;
};

template < template < class, class, int > class AlgorithmType, class DataType, int DIM >
class StandardBoundaryScalarAssembly
    : public AssemblyAlgorithmBase<AlgorithmType, StandardBoundaryScalarAssembly< AlgorithmType, DataType, DIM >, DataType, DIM > 
{
public:
  typedef std::vector< DataType > LocalObjectType;
  typedef hiflow::Quadrature< DataType > QuadratureType;
  typedef AssemblyAlgorithmBase< AlgorithmType, StandardBoundaryScalarAssembly< AlgorithmType, DataType, DIM >, DataType, DIM > Base;

  // Name resolution does not manage to get the base members, so
  // we must do it ourselves.
  using Base::curr_;
  using Base::elem_;
  using Base::has_next;
  using Base::space_;
  using Base::traversal_;

  StandardBoundaryScalarAssembly(const VectorSpace< DataType, DIM > &space,
                                 std::vector< DataType > &values)
      : Base(space), values_(values) {
    remove_non_local_elements();
    sort_elements(space, this->traversal_);
  }

  /// The has_next() and next() functions are overloaded in
  /// order to skip elements that do not belong to the local
  /// subdomain. This is done by setting the corresponding
  /// entries in the traversal_ array to -1, and later skipping
  /// those items.

  const Element< DataType, DIM > &next() {
    assert(this->has_next());

    this->elem_ =
        Element< DataType, DIM >(this->space_, this->traversal_[this->curr_]);

    ++(this->curr_);

    return this->elem_;
  }

  void add(const Element< DataType, DIM > &element,
           const LocalObjectType &local_val) {
    mesh::TDim tdim = this->space_.mesh().tdim();
    mesh::IncidentEntityIterator iter =
        element.get_cell().begin_incident(tdim - 1);
    mesh::IncidentEntityIterator end =
        element.get_cell().end_incident(tdim - 1);
    int facet_number = 0;
    for (; iter != end; iter++) {

      this->values_[iter->id()] += local_val[facet_number];
      ++facet_number;
    }
  }

  void reset(typename GlobalAssembler< DataType, DIM >::LocalVector &local_vec) {
    mesh::TDim tdim = this->space_.mesh().tdim();
    local_vec.clear();
    local_vec.resize(this->elem_.get_cell().num_incident_entities(tdim - 1),
                     0.);
  }

private:
  /// Remove non_local elements

  void remove_non_local_elements() {

    const mesh::Mesh &mesh = this->space_.mesh();

    if (!mesh.has_attribute("_remote_index_", mesh.tdim())) {
      // If the "_remote_index_" attribute does not exist, we
      // assume that there are no ghost cells.
      return;
    }

    int remote_index;
    int index;
    std::vector< int > traversal_tmp;
    traversal_tmp.reserve(mesh.num_entities(mesh.tdim()));

    for (mesh::EntityIterator it_cell = mesh.begin(mesh.tdim()),
                              e_it_cell = mesh.end(mesh.tdim());
         it_cell != e_it_cell; ++it_cell) {
      // test if cell on subdomain
      it_cell->get("_remote_index_", &remote_index);
      if (remote_index == -1) {
        index = it_cell->index();
        traversal_tmp.push_back(index);
      }
    }
    this->traversal_ = traversal_tmp;
  }

  std::vector< DataType > &values_;
};


template < template < class, class, int > class AlgorithmType, class DataType, int DIM >
class StandardBoundaryMultipleScalarAssembly
    : public AssemblyAlgorithmBase<
          AlgorithmType,
          StandardBoundaryMultipleScalarAssembly< AlgorithmType, DataType, DIM >,
          DataType,
          DIM > {
public:
  typedef std::vector< std::vector<DataType> > LocalObjectType;
  typedef hiflow::Quadrature< DataType > QuadratureType;
  typedef AssemblyAlgorithmBase<
      AlgorithmType, StandardBoundaryMultipleScalarAssembly< AlgorithmType, DataType, DIM >,
      DataType,
      DIM >
      Base;

  // Name resolution does not manage to get the base members, so
  // we must do it ourselves.
  using Base::curr_;
  using Base::elem_;
  using Base::has_next;
  using Base::space_;
  using Base::traversal_;

  StandardBoundaryMultipleScalarAssembly(const VectorSpace< DataType, DIM > &space,
                                         std::vector< std::vector<DataType> > &values,
                                         size_t num_scalars)
      : Base(space), values_(values), num_scalars_(num_scalars) 
  {
    this->remove_non_local_elements();
    sort_elements(space, this->traversal_);
  }

  /// The has_next() and next() functions are overloaded in
  /// order to skip elements that do not belong to the local
  /// subdomain. This is done by setting the corresponding
  /// entries in the traversal_ array to -1, and later skipping
  /// those items.

  const Element< DataType, DIM > &next() {
    assert(this->has_next());

    this->elem_ = Element< DataType, DIM >(this->space_, this->traversal_[this->curr_]);

    ++(this->curr_);

    return this->elem_;
  }

  void add(const Element< DataType, DIM > &element,
           const LocalObjectType &local_val) 
  {
    mesh::TDim tdim = this->space_.mesh().tdim();
    mesh::IncidentEntityIterator iter = element.get_cell().begin_incident(tdim - 1);
    mesh::IncidentEntityIterator end = element.get_cell().end_incident(tdim - 1);
    int facet_number = 0;
    for (; iter != end; iter++) 
    {
      const int facet_id = iter->id();
      assert (local_val[facet_number].size() == this->num_scalars_);
      assert (this->values_[facet_id].size() == this->num_scalars_);
      
      for (size_t l=0; l<this->num_scalars_; ++l)
      {
        this->values_[facet_id][l] += local_val[facet_number][l];
      }
      ++facet_number;
    }
  }

  void reset(std::vector<typename GlobalAssembler< DataType, DIM >::LocalVector> &local_vec) 
  {
    mesh::TDim tdim = this->space_.mesh().tdim();
    local_vec.clear();
    local_vec.resize(this->elem_.get_cell().num_incident_entities(tdim - 1));
    for (int f=0; f<local_vec.size(); ++f)
    {
      local_vec[f].resize(this->num_scalars_, 0.);
    }
  }

private:
  /// Remove non_local elements

  void remove_non_local_elements() {

    const mesh::Mesh &mesh = this->space_.mesh();

    if (!mesh.has_attribute("_remote_index_", mesh.tdim())) {
      // If the "_remote_index_" attribute does not exist, we
      // assume that there are no ghost cells.
      return;
    }

    int remote_index;
    int index;
    std::vector< int > traversal_tmp;
    traversal_tmp.reserve(mesh.num_entities(mesh.tdim()));

    for (mesh::EntityIterator it_cell = mesh.begin(mesh.tdim()),
                              e_it_cell = mesh.end(mesh.tdim());
         it_cell != e_it_cell; ++it_cell) {
      // test if cell on subdomain
      it_cell->get("_remote_index_", &remote_index);
      if (remote_index == -1) {
        index = it_cell->index();
        traversal_tmp.push_back(index);
      }
    }
    this->traversal_ = traversal_tmp;
  }
  
  std::vector< std::vector<DataType> > &values_;
  size_t num_scalars_;
};

//////////////// end helper functions ////////////////

//////////////// Implementation of StandardGlobalAssembler ////////////////
template < class DataType, int DIM >
class StandardGlobalAssembler : public GlobalAssembler< DataType, DIM > 
{
  typedef GlobalAssembler< DataType, DIM > GlobalAsm;
  typedef VectorSpace< DataType, DIM > VecSpace;

  virtual void assemble_scalar_impl(const VecSpace &space,
                                    typename GlobalAsm::ScalarAssemblyFunction local_asm,
                                    std::vector< DataType > &vec,
                                    typename GlobalAsm::QuadratureSelectionFunction q_select) const;

  virtual void assemble_multiple_scalar_impl(const VecSpace &space,
                                             typename GlobalAsm::MultipleScalarAssemblyFunction local_asm,
                                             int num_scalars, 
                                             std::vector< std::vector< DataType > > &vec,
                                             typename GlobalAsm::QuadratureSelectionFunction q_select) const;

  virtual void assemble_vector_impl(const VecSpace &space,
                                    typename GlobalAsm::VectorAssemblyFunction local_asm,
                                    typename GlobalAsm::GlobalVector &vec,
                                    typename GlobalAsm::QuadratureSelectionFunction q_select) const;

  virtual void assemble_matrix_impl(const VecSpace &space,
                                    typename GlobalAsm::MatrixAssemblyFunction local_asm,
                                    typename GlobalAsm::GlobalMatrix &mat,
                                    typename GlobalAsm::QuadratureSelectionFunction q_select) const;

  virtual void assemble_scalar_boundary_impl(const VecSpace &space,
                                             typename GlobalAsm::BoundaryScalarAssemblyFunction local_asm,
                                             std::vector< DataType > &vec,
                                             typename GlobalAsm::FacetQuadratureSelectionFunction fq_select) const;

  virtual void assemble_multiple_scalar_boundary_impl(const VecSpace &space,
                                                      typename GlobalAsm::BoundaryMultipleScalarAssemblyFunction local_asm,
                                                      int num_scalars, 
                                                      std::vector< std::vector< DataType > > &vec,
                                                      typename GlobalAsm::FacetQuadratureSelectionFunction fq_select) const;
                                                     
  virtual void assemble_vector_boundary_impl(const VecSpace &space,
                                             typename GlobalAsm::BoundaryVectorAssemblyFunction local_asm,
                                             typename GlobalAsm::GlobalVector &vec,
                                             typename GlobalAsm::FacetQuadratureSelectionFunction fq_select) const;

  virtual void assemble_matrix_boundary_impl(const VecSpace &space,
                                             typename GlobalAsm::BoundaryMatrixAssemblyFunction local_asm,
                                             typename GlobalAsm::GlobalMatrix &mat,
                                             typename GlobalAsm::FacetQuadratureSelectionFunction fq_select) const;
};

template < class DataType, int DIM >
void StandardGlobalAssembler< DataType, DIM >::assemble_scalar_impl(const VecSpace &space,
                                                                    typename GlobalAsm::ScalarAssemblyFunction local_asm,
                                                                    std::vector< DataType > &vec,
                                                                    typename GlobalAsm::QuadratureSelectionFunction q_select) const 
{
  StandardScalarAssembly< InteriorAssemblyAlgorithm, DataType, DIM > assembly(space, vec);
  assembly.assemble(local_asm, q_select);
}

template < class DataType, int DIM >
void StandardGlobalAssembler< DataType, DIM >::assemble_multiple_scalar_impl(const VecSpace &space,
                                                                             typename GlobalAsm::MultipleScalarAssemblyFunction local_asm,
                                                                             const int num_scalars, std::vector< std::vector< DataType > > &vec,
                                                                             typename GlobalAsm::QuadratureSelectionFunction q_select) const 
{
  StandardMultipleScalarAssembly< InteriorAssemblyAlgorithm, DataType, DIM > assembly(space, vec, num_scalars);
  assembly.assemble(local_asm, q_select);
}

template < class DataType, int DIM >
void StandardGlobalAssembler< DataType, DIM >::assemble_vector_impl(const VecSpace &space,
                                                                    typename GlobalAsm::VectorAssemblyFunction local_asm,
                                                                    typename GlobalAsm::GlobalVector &vec,
                                                                    typename GlobalAsm::QuadratureSelectionFunction q_select) const 
{
  StandardVectorAssembly< InteriorAssemblyAlgorithm, DataType, DIM > assembly(space, vec);
  assembly.assemble(local_asm, q_select);
}

template < class DataType, int DIM >
void StandardGlobalAssembler< DataType, DIM >::assemble_matrix_impl(const VecSpace &space,
                                                                    typename GlobalAsm::MatrixAssemblyFunction local_asm,
                                                                    typename GlobalAsm::GlobalMatrix &mat,
                                                                    typename GlobalAsm::QuadratureSelectionFunction q_select) const 
{
  StandardMatrixAssembly< InteriorAssemblyAlgorithm, DataType, DIM > assembly(space, mat);
  assembly.assemble(local_asm, q_select);
}

template < class DataType, int DIM >
void StandardGlobalAssembler< DataType, DIM >::assemble_scalar_boundary_impl(const VecSpace &space,
                                                                             typename GlobalAsm::BoundaryScalarAssemblyFunction local_asm,
                                                                             std::vector< DataType > &vec,
                                                                             typename GlobalAsm::FacetQuadratureSelectionFunction fq_select) const 
{
  // TODO: how should vec be defined -> all facets or only boundary facets?
  // If the latter, then what ordering should be used?
  StandardBoundaryScalarAssembly< BoundaryAssemblyAlgorithm, DataType, DIM > assembly(space, vec);
  assembly.assemble(local_asm, fq_select);
}

template < class DataType, int DIM >
void StandardGlobalAssembler< DataType, DIM >::assemble_multiple_scalar_boundary_impl(const VecSpace &space,
                                                                                      typename GlobalAsm::BoundaryMultipleScalarAssemblyFunction local_asm,
                                                                                      const int num_scalars,
                                                                                      std::vector< std::vector<DataType> > &vec,
                                                                                      typename GlobalAsm::FacetQuadratureSelectionFunction fq_select) const 
{
  // TODO: how should vec be defined -> all facets or only boundary facets?
  // If the latter, then what ordering should be used?
  StandardBoundaryMultipleScalarAssembly< BoundaryAssemblyAlgorithm, DataType, DIM > assembly(space, vec, num_scalars);
  assembly.assemble(local_asm, fq_select);
}

template < class DataType, int DIM >
void StandardGlobalAssembler< DataType, DIM >::assemble_vector_boundary_impl(const VecSpace &space,
                                                                             typename GlobalAsm::BoundaryVectorAssemblyFunction local_asm,
                                                                             typename GlobalAsm::GlobalVector &vec,
                                                                             typename GlobalAsm::FacetQuadratureSelectionFunction fq_select) const 
{
  StandardVectorAssembly< BoundaryAssemblyAlgorithm, DataType, DIM > assembly(space, vec);
  assembly.assemble(local_asm, fq_select);
}

template < class DataType, int DIM >
void StandardGlobalAssembler< DataType, DIM >::assemble_matrix_boundary_impl(const VecSpace &space,
                                                                             typename GlobalAsm::BoundaryMatrixAssemblyFunction local_asm,
                                                                             typename GlobalAsm::GlobalMatrix &mat,
                                                                             typename GlobalAsm::FacetQuadratureSelectionFunction fq_select) const 
{
  StandardMatrixAssembly< BoundaryAssemblyAlgorithm, DataType, DIM > assembly(space, mat);
  assembly.assemble(local_asm, fq_select);
}

} // namespace hiflow
#endif
