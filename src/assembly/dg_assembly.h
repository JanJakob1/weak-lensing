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

#ifndef HIFLOW_ASSEMBLY_DG_ASSEMBLY_H
#define HIFLOW_ASSEMBLY_DG_ASSEMBLY_H

#include "assembly/global_assembler.h"
#include "assembly/standard_assembly.h"
#include "assembly/quadrature_selection.h"
#include "mesh/interface.h"
#include "quadrature/quadrature.h"

#include "common/log.h"
#include "common/sort_permutation.h"
#include "common/sorted_array.h"
#include "mesh/mesh.h"
#include <boost/function.hpp>
#include <functional>

// this seems to be necessary, at least when using Elements based on facet moments (e.g. BDM, RT)
#define INCLUDE_GHOST_INTERFACES

/// \author Staffan Ronnas, Jonathan Schwegler, Simon Gawlok
const int PRINT_RANK = 0;

namespace hiflow {

/// \author Staffan Ronnas, Jonathan Schwegler, Simon Gawlok

/// \brief Class for performing global assembly over interfaces for e.g.
/// Discontinuous Galerkin-Methods.

template < class DataType, int DIM >
class DGGlobalAssembler : public StandardGlobalAssembler< DataType, DIM > {
public:
  enum InterfaceSide {
    INTERFACE_MASTER = 0,
    INTERFACE_SLAVE = 1,
    INTERFACE_BOUNDARY = 2,
    INTERFACE_NONE = -1
  };
  
  // TODO: this could be optimized by splitting into separate cell/facet
  // selection functions.
  typedef std::function < void ( const Element< DataType, DIM > &, const Element< DataType, DIM > &, 
                                 int, int,
                                 Quadrature< DataType > &, Quadrature< DataType > &) >
                                 IFQuadratureSelectionFun;
  
  typedef GlobalAssembler< DataType, DIM > GlobalAsm;
  typedef VectorSpace< DataType, DIM > VecSpace;

  DGGlobalAssembler();

  ///
  /// \brief Assemble global matrix for a VecSpace.
  ///
  /// \details Assembles a matrix
  /// \f$A_{ij} = \sum_{K \in \mathcal T_h}\int_{\partial K}{f(x,\varphi_i,
  /// \varphi_j)dx}\f$ over the interfaces defined by the mesh associated to a
  /// VecSpace. The integrand is defined through a
  /// InterfaceMatrixAssemblyFun, which should return the locally assembled
  /// matrix for each element.
  ///
  /// \param[in] space           the VecSpace for which the assembly is
  /// performed \param[in] local_asm       function or functor that performs
  /// local matrix assembly \param[out] global_matrix  the assembled matrix
  /// \f$A_{ij}\f$ \see concept_assembly
  ///
  /// template Argument LocalAssembler should provide the following routine:
  /// void operator() (const Element< DataType, DIM > &master_elem,
  ///                  const Element< DataType, DIM > &slave_elem,
  ///                  const Quadrature< DataType > &master_quad,
  ///                  const Quadrature< DataType > &slave_quad,
  ///                  int master_facet_number, 
  ///                  int slave_facet_number,
  ///                  InterfaceSide if_side, 
  ///                  int slave_index, int num_slaves, 
  ///                  typename GlobalAssembler< DataType, DIM >::LocalMatrix &vals)
  
  template<class LocalAssembler>
  void assemble_interface_matrix(const VecSpace &space,
                                 LocalAssembler& local_asm,
                                 typename GlobalAsm::GlobalMatrix &matrix) const;

  ///
  /// \brief Assemble global vector for a VecSpace.
  ///
  /// \details Assembles a vector
  /// \f$b_i = \sum_{K \in \mathcal T_h} \int_{\partial K}{f(x, \varphi_i)dx}\f$
  /// over the interfaces defined by the mesh
  /// associated to a VecSpace. The integrand is defined through
  /// a InterfaceVectorAssemblyFun, which should return the locally
  /// assembled vector for each element.
  ///
  /// \param[in] space      the VecSpace for which the assembly is performed
  /// \param[in] local_asm  function or functor that performs local vector
  /// assembly \param[out] global_vector  the assembled vector \f$b_i\f$ \see
  /// concept_assembly
  ///
  /// template Argument LocalAssembler should provide the following routine:
  /// void operator() (const Element< DataType, DIM > &master_elem,
  ///                  const Element< DataType, DIM > &slave_elem,
  ///                  const Quadrature< DataType > &master_quad,
  ///                  const Quadrature< DataType > &slave_quad,
  ///                  int master_facet_number, 
  ///                  int slave_facet_number,
  ///                  InterfaceSide if_side, 
  ///                  int slave_index, int num_slaves, 
  ///                  typename GlobalAssembler< DataType, DIM >::LocalVector &vals)
  
  template<class LocalAssembler>
  void assemble_interface_vector(const VecSpace &space,
                                 LocalAssembler& local_asm,
                                 typename GlobalAsm::GlobalVector &vec) const;

  ///
  /// \brief Assemble global value for a VecSpace.
  ///
  /// \details Assembles a scalar
  /// \f$v_K = \int_{\partial K}{f(x, \varphi_i)dx}, \ K \in \mathcal T_h\f$
  /// over the interfaces defined by the mesh
  /// associated to a VecSpace. The integrand is defined through
  /// a InterfaceScalarAssemblyFun, which should return the locally
  /// assembled value for each element.
  ///
  /// \param[in] space      the VecSpace for which the assembly is performed
  /// \param[in] local_asm  function or functor that performs local vector
  /// assembly \param[out] values  the assembled values \f$v_K\f$ \see
  /// concept_assembly
  ///
  /// template Argument LocalAssembler should provide the following routine:
  /// void operator() (const Element< DataType, DIM > &master_elem,
  ///                  const Element< DataType, DIM > &slave_elem,
  ///                  const Quadrature< DataType > &master_quad,
  ///                  const Quadrature< DataType > &slave_quad,
  ///                  int master_facet_number, 
  ///                  int slave_facet_number,
  ///                  InterfaceSide if_side, 
  ///                  int slave_index, int num_slaves, 
  ///                  DataType &value)
                                   
  template<class LocalAssembler>
  void assemble_interface_scalar(const VecSpace &space,
                                 LocalAssembler& local_asm,
                                 std::vector< DataType > &values) const;

  ///
  /// \brief Assemble global value for a VecSpace.
  ///
  /// \details Assembles a scalar
  /// \f_K = \int_{\partial K}{f(x)dx}, \ K \in \mathcal T_h\f$
  /// over the interfaces defined by the mesh
  /// associated to a VecSpace. The integrand is defined through
  /// a InterfaceScalarAssemblyFun, which should return the locally
  /// assembled value for each interface.
  ///
  /// \param[in] space      the VecSpace for which the assembly is performed
  /// \param[in] local_asm  function or functor that performs local vector
  /// assembly \param[out] values  the assembled values \f$v_K\f$ correctly
  /// distributed to cells \see concept_assembly
  ///
  ///
  /// template Argument LocalAssembler should provide the following routine:
  /// void operator() (const Element< DataType, DIM > &master_elem,
  ///                  const Element< DataType, DIM > &slave_elem,
  ///                  const Quadrature< DataType > &master_quad,
  ///                  const Quadrature< DataType > &slave_quad,
  ///                  int master_facet_number, 
  ///                  int slave_facet_number,
  ///                  InterfaceSide if_side, 
  ///                  int slave_index, int num_slaves, 
  ///                  DataType &value)
                                   
  template<class LocalAssembler>
  void assemble_interface_scalar_cells(const VecSpace &space,
                                       LocalAssembler& local_asm,
                                       std::vector< DataType > &values) const;

  ///
  /// \brief Assemble global value for a VecSpace.
  ///
  /// \details Assembles multiple scalars
  /// \f_K = \int_{\partial K}{f(x)dx}, \ K \in \mathcal T_h\f$
  /// over the interfaces defined by the mesh
  /// associated to a VecSpace. The integrand is defined through
  /// a InterfaceMultipleScalarAssemblyFun, which should return the locally
  /// assembled values for each interface.
  ///
  /// \param[in] space      the VecSpace for which the assembly is performed
  /// \param[in] local_asm  function or functor that performs local vector
  /// assembly \param[in] num_scalars dimension of integrand \param[out] values
  /// the assembled values \f_K\f$ correctly distributed to cells \see
  /// concept_assembly
  ///
  /// template Argument LocalAssembler should provide the following routine:
  /// void operator() (const Element< DataType, DIM > &master_elem,
  ///                  const Element< DataType, DIM > &slave_elem,
  ///                  const Quadrature< DataType > &master_quad,
  ///                  const Quadrature< DataType > &slave_quad,
  ///                  int master_facet_number, 
  ///                  int slave_facet_number,
  ///                  InterfaceSide if_side, 
  ///                  int slave_index, int num_slaves, 
  ///                  std::vector<DataType> &vals)
                                   
  template<class LocalAssembler>
  void assemble_interface_multiple_scalar_cells(const VecSpace &space,
                                                LocalAssembler& local_asm, 
                                                const int num_scalars,
                                                std::vector< typename GlobalAsm::LocalVector > &values) const;

  ///
  /// \brief Add contributions of an interface scalar assembly to a vector of
  /// cell values in a naive way. This functionality is, e.g., needed in error
  /// estimators
  void distribute_interface_to_cell_values_naive(
      const VecSpace &space,
      std::vector< DataType > &cell_values,
      const std::vector< DataType > &interface_values) const;

  ///
  /// \brief Set the function used to determine which quadrature rule should be
  /// used. for interface quadrature.
  ///
  /// \details The choice of quadrature rules can be controlled
  /// by providing a IFQuadratureSelectionFunction. This function
  /// or functor will be called before the local assembly is
  /// performed on each element. By default, the
  /// DefaultInterfaceQuadratureSelection function is used.
  ///
  /// \param[in] q_select   new quadrature selection function
  void
  set_interface_quadrature_selection_fun(IFQuadratureSelectionFun q_select);

private:


  IFQuadratureSelectionFun if_q_select_;
};

// Add entries in local matrix to global matrix.

template < class DataType, int DIM >
void add_global(const VectorSpace< DataType, DIM > &space,
                const int test_cell_index,
                const int trial_cell_index,
                const std::vector< int > &row_dofs,
                const std::vector< int > &col_dofs,
                const typename GlobalAssembler< DataType, DIM >::LocalMatrix &lm,
                typename GlobalAssembler< DataType, DIM >::GlobalMatrix &gm) {

  const size_t num_dofs_row = row_dofs.size();
  const size_t num_dofs_col = col_dofs.size();

  std::vector< int > row_dofs_sort_permutation;
  std::vector< int > col_dofs_sort_permutation;
  std::vector< int > row_dofs_sorted(num_dofs_row);
  std::vector< int > col_dofs_sorted(num_dofs_col);

  // get permutation for sorting dofs
  sortingPermutation(row_dofs, row_dofs_sort_permutation);
  sortingPermutation(col_dofs, col_dofs_sort_permutation);

  // fill sorted dof array
  for (size_t i = 0; i != num_dofs_row; ++i) {
    row_dofs_sorted[i] = row_dofs[row_dofs_sort_permutation[i]];
  }
  for (size_t i = 0; i != num_dofs_col; ++i) {
    col_dofs_sorted[i] = col_dofs[col_dofs_sort_permutation[i]];
  }

  std::vector< DataType > dof_factors_trial;
  std::vector< DataType > dof_factors_test;
  
  space.dof().get_dof_factors_on_cell (trial_cell_index, dof_factors_trial);
  space.dof().get_dof_factors_on_cell (test_cell_index, dof_factors_test);
  
  assert (dof_factors_trial.size() == num_dofs_col);
  assert (dof_factors_test.size() == num_dofs_row);
    
  // create row array
  std::vector< int > row_indices;
  std::vector< int > row_permutation;
  row_indices.reserve(num_dofs_row);
  row_permutation.reserve(num_dofs_row);
  for (size_t i = 0; i != num_dofs_row; ++i) {
    if (space.dof().is_dof_on_subdom(row_dofs_sorted[i])) {
      row_indices.push_back(row_dofs_sorted[i]);
      row_permutation.push_back(row_dofs_sort_permutation[i]);
    }
  }

  // fill reduced and sorted local matrix
  typename GlobalAssembler< DataType, DIM >::LocalMatrix local_mat_sorted_reduced;
  if (!row_indices.empty() && col_dofs_sorted.size() > 0) 
  {
    local_mat_sorted_reduced.Resize(row_indices.size(), col_dofs_sorted.size());
    for (size_t i = 0; i != row_indices.size(); ++i) 
    {
      const int test_ind = row_permutation[i];
      for (size_t j = 0; j != col_dofs_sorted.size(); ++j) 
      {
        const int trial_ind = col_dofs_sort_permutation[j];
        local_mat_sorted_reduced(i, j) = dof_factors_trial[trial_ind]
                                       * dof_factors_test[test_ind]
                                       * lm(test_ind, trial_ind);
      }
    }

    // Add local to global matrix
    gm.Add(vec2ptr(row_indices), row_indices.size(), vec2ptr(col_dofs_sorted),
           col_dofs_sorted.size(), &local_mat_sorted_reduced(0, 0));
  }
}

template < class DataType, int DIM >
void add_global(const VectorSpace< DataType, DIM > &space,
                const int test_cell_index,
                const std::vector< int > &dofs,
                const typename GlobalAssembler< DataType, DIM >::LocalVector &lv,
                typename GlobalAssembler< DataType, DIM >::GlobalVector &vec) {

  const size_t num_dofs = dofs.size();

  std::vector< int > dofs_sort_permutation;
  std::vector< int > dofs_sorted(num_dofs);

  // get permutation for sorting dofs
  sortingPermutation(dofs, dofs_sort_permutation);

  // fill sorted dof array
  for (size_t i = 0; i != num_dofs; ++i) {
    dofs_sorted[i] = dofs[dofs_sort_permutation[i]];
  }

  std::vector< DataType > dof_factors_test;
  space.dof().get_dof_factors_on_cell (test_cell_index, dof_factors_test);
  
  assert (dof_factors_test.size() == num_dofs);


  // create row array
  std::vector< int > row_indices;
  row_indices.reserve(num_dofs);

  typename GlobalAssembler< DataType, DIM >::LocalVector local_vec_sorted;
  local_vec_sorted.reserve(num_dofs);

  for (size_t i = 0; i != num_dofs; ++i) 
  {
    if (space.dof().is_dof_on_subdom(dofs_sorted[i])) 
    {
      const int test_ind = dofs_sort_permutation[i];
      row_indices.push_back(dofs_sorted[i]);
      local_vec_sorted.push_back(dof_factors_test[test_ind] * lv[test_ind]);
    }
  }

  // Add local to global vector
  if (!row_indices.empty()) {
    vec.Add(vec2ptr(row_indices), row_indices.size(),
            vec2ptr(local_vec_sorted));
  }
}

template < class DataType, int DIM >
DGGlobalAssembler< DataType, DIM >::DGGlobalAssembler()
    : if_q_select_(DefaultInterfaceQuadratureSelection< DataType, DIM >()) 
{
  // By default, we don't reset the target.
  this->should_reset_assembly_target(false);
}

template < class DataType, int DIM >
template< class LocalAssembler>
void DGGlobalAssembler< DataType, DIM >::assemble_interface_matrix( const VecSpace &space, 
                                                                    LocalAssembler& local_asm,
                                                                    typename GlobalAsm::GlobalMatrix &matrix) const 
{
  if (this->should_reset_assembly_target_) 
  {
    matrix.Zeros();
  }
  
  // Create interface list from mesh
  mesh::ConstMeshPtr mesh = &space.mesh(); // OK since pointer is intrusive
  mesh::InterfaceList if_list = mesh::InterfaceList::create(mesh);

  typename GlobalAsm::LocalMatrix L_MM, L_MS, L_SM, L_SS;

  std::vector< int > master_dofs, slave_dofs;

  Quadrature< DataType > master_quadrature, slave_quadrature;

  // Loop over interfaces

  for (mesh::InterfaceList::const_iterator it = if_list.begin(),
       end_it = if_list.end();
       it != end_it; ++it) 
  {
    int remote_index_master = -10;
    mesh->get_attribute_value("_remote_index_", mesh->tdim(),
                              it->master_index(), &remote_index_master);

    // Master dofs
    const int master_cell_index = it->master_index();
    Element< DataType, DIM > master_elem(space, master_cell_index);
    
    space.get_dof_indices(master_cell_index, &master_dofs);
    const int master_facet_number = it->master_facet_number();

    L_MM.Resize(master_dofs.size(), master_dofs.size());
    L_MM.Zeros();

    // Initialize master quadrature
    if_q_select_(master_elem, master_elem, master_facet_number,
                 master_facet_number, master_quadrature, master_quadrature);

    const int num_slaves = it->num_slaves(); 
    assert (num_slaves >= 0);
   
   // subdom              = localdom + ghostdom
   // num_slaves = 0      -> bdy of subdom
   // num_slaves > 0      -> interior of subdom
   // remote_index == -1  -> cell \in localdom
   // remote_index >= 0   -> cell \in ghostdom
   
   // num_slaves = 0 / remote_index_master == -1 -> bdy of localdom = physical bdy
   // num_slaves = 0 / remote_index_master >= 0  -> bdy of ghost    = physical interior
   
    // treat boundary facet
    if (remote_index_master == -1) 
    {
      if (num_slaves == 0) 
      {
        local_asm(master_elem, master_elem, 
                  master_quadrature, master_quadrature, 
                  master_facet_number, master_facet_number,
                  INTERFACE_MASTER, INTERFACE_BOUNDARY, 
                  -1, num_slaves,
                  L_MM);
        assert (!L_MM.contains_nan());
      }
      add_global(space, master_cell_index, master_cell_index, master_dofs, master_dofs, L_MM, matrix);
    }
    
    // treat interior facets
    // Loop over slaves
    for (int s = 0; s < num_slaves; ++s) 
    {
      const int slave_cell_index = it->slave_index(s);
      const int slave_facet_number = it->slave_facet_number(s);

      int remote_index_slave = -10;
      mesh->get_attribute_value("_remote_index_", mesh->tdim(),
                                slave_cell_index, &remote_index_slave);
                                
      Element< DataType, DIM > slave_elem(space, slave_cell_index);
      space.get_dof_indices(slave_cell_index, &slave_dofs);
            

      // Initialize slave quadrature. NB: only called once per slave.
      // default quad selection: all quadrature points lie on that part of the interface,
      // which has non-empty intersection with slave cell
      if_q_select_(master_elem, slave_elem, 
                   master_facet_number, slave_facet_number, 
                   master_quadrature, slave_quadrature);

#ifdef INCLUDE_GHOST_INTERFACES
      bool add_to_master = true;
      bool add_to_slave = true;
#else
      bool add_to_master = (remote_index_master == -1); 
      bool add_to_slave = (remote_index_slave == -1); 
#endif

      if (add_to_master) 
      {
        // master / slave
        L_MS.Resize(master_dofs.size(), slave_dofs.size());
        L_MS.Zeros();
        local_asm(master_elem, slave_elem, 
                  master_quadrature, slave_quadrature,
                  master_facet_number, slave_facet_number, 
                  INTERFACE_SLAVE, INTERFACE_MASTER, 
                  s, num_slaves,
                  L_MS);
        assert (!L_MS.contains_nan());
        add_global(space, master_cell_index, slave_cell_index, master_dofs, slave_dofs, L_MS, matrix);
      
        // master / master
        // Note: in case of hanging nodes, there holds num_slaves > 1, 
        // i.e. the combination master / master is called more then once.
        // Use the slave index s passed to your local assembler implementation 
        // to decide what to do in this case.
        // Why do we do that? Because for evaluating facet integrals, it might be necessary to
        // have both master and slave element at hand, even if the combination master / master is considered. 
        // Example: need to compute jump and average of discontinuous function, which is neither ansatz, nor 
        // test function
        //
        // Note: actually, calling master / master more than once should be fine,
        // since the quadrature is restricted to the slave portion of the interface
        
        L_MM.Resize(master_dofs.size(), master_dofs.size());
        L_MM.Zeros();
        local_asm(master_elem, slave_elem, 
                  master_quadrature, slave_quadrature,
                  master_facet_number, slave_facet_number, 
                  INTERFACE_MASTER, INTERFACE_MASTER, 
                  s, num_slaves,
                  L_MM);
        assert (!L_MM.contains_nan());
        add_global(space, master_cell_index, master_cell_index, master_dofs, master_dofs, L_MM, matrix);
      }
      if (add_to_slave) 
      {
        // slave / master
        L_SM.Resize(slave_dofs.size(), master_dofs.size());
        L_SM.Zeros();
        local_asm(master_elem, slave_elem,  
                  master_quadrature, slave_quadrature,
                  master_facet_number, slave_facet_number, 
                  INTERFACE_MASTER, INTERFACE_SLAVE, 
                  s, num_slaves,
                  L_SM);
        assert (!L_SM.contains_nan());
        add_global(space, slave_cell_index, master_cell_index, slave_dofs, master_dofs, L_SM, matrix);

        // slave / slave
        L_SS.Resize(slave_dofs.size(), slave_dofs.size());
        L_SS.Zeros();
        local_asm(master_elem, slave_elem, 
                  master_quadrature, slave_quadrature,
                  master_facet_number, slave_facet_number, 
                  INTERFACE_SLAVE, INTERFACE_SLAVE, 
                  s, num_slaves,
                  L_SS);
        assert (!L_SS.contains_nan());
        add_global(space, slave_cell_index, slave_cell_index, slave_dofs, slave_dofs, L_SS, matrix);
      }
    }
  }
}

template < class DataType, int DIM >
template< class LocalAssembler>
void DGGlobalAssembler< DataType, DIM >::assemble_interface_vector(const VecSpace &space, 
                                                                   LocalAssembler& local_asm, 
                                                                   typename GlobalAsm::GlobalVector &vec) const 
{
  if (this->should_reset_assembly_target_) 
  {
    vec.Zeros();
  }
  // Create interface list from mesh
  mesh::ConstMeshPtr mesh = &space.mesh(); // OK since pointer is intrusive

  mesh::InterfaceList if_list = mesh::InterfaceList::create(mesh);

  typename GlobalAsm::LocalVector L_M, L_S;

  std::vector< int > master_dofs, slave_dofs;

  Quadrature< DataType > master_quadrature, slave_quadrature;

  // Loop over interfaces

  for (mesh::InterfaceList::const_iterator it = if_list.begin(),
       end_it = if_list.end();
       it != end_it; ++it) 
  {
    const int master_cell_index = it->master_index();
    const int master_facet_number = it->master_facet_number();
    int remote_index_master = -10;
    mesh->get_attribute_value("_remote_index_", mesh->tdim(),
                              master_cell_index, &remote_index_master);

    // Master dofs
    Element< DataType, DIM > master_elem(space, master_cell_index);
    space.get_dof_indices(master_cell_index, &master_dofs);
    
    L_M.clear();
    L_M.resize(master_dofs.size(), 0.);

    // Initialize master quadrature
    if_q_select_(master_elem, master_elem, 
                 master_facet_number, master_facet_number, 
                 master_quadrature, master_quadrature);

    const int num_slaves = it->num_slaves();
    if (remote_index_master == -1) 
    {
      if (num_slaves == 0) 
      {
        // boundary facet
        local_asm(master_elem, master_elem, 
                  master_quadrature, master_quadrature, 
                  master_facet_number, master_facet_number,
                  INTERFACE_BOUNDARY, 
                  -1, num_slaves,
                  L_M);
        assert (!contains_nan(L_M));
      }
      add_global(space, master_cell_index, master_dofs, L_M, vec);
    }
    // Loop over slaves
    for (int s = 0; s < num_slaves; ++s) 
    {
      const int slave_facet_number = it->slave_facet_number(s);
      const int slave_cell_index = it->slave_index(s);
      
      int remote_index_slave = -10;
      mesh->get_attribute_value("_remote_index_", mesh->tdim(),
                                slave_cell_index, &remote_index_slave);
                                
      Element< DataType, DIM > slave_elem(space, slave_cell_index);
      space.get_dof_indices(slave_cell_index, &slave_dofs);

      
      // Initialize slave quadrature. NB: only called once per slave.
      // default quad selection: all quadrature points lie on that part of the interface,
      // which has non-empty intersection with slave cell
      if_q_select_(master_elem, slave_elem, 
                   master_facet_number, slave_facet_number, 
                   master_quadrature, slave_quadrature);

#ifdef INCLUDE_GHOST_INTERFACES
      bool add_to_master = true;
      bool add_to_slave = true;
#else
      bool add_to_master = (remote_index_master == -1); 
      bool add_to_slave = (remote_index_slave == -1); 
#endif

      if (add_to_master) 
      {
        // master
        L_M.clear();
        L_M.resize(master_dofs.size(), 0.);
        local_asm(master_elem, slave_elem, 
                  master_quadrature, slave_quadrature, 
                  master_facet_number, slave_facet_number,
                  INTERFACE_MASTER, 
                  s, num_slaves,
                  L_M);
        assert (!contains_nan(L_M));
        add_global(space, master_cell_index, master_dofs, L_M, vec);
      }
      if (add_to_slave) 
      {
        // slave
        L_S.clear();
        L_S.resize(slave_dofs.size(), 0.);
        local_asm(master_elem, slave_elem, 
                  master_quadrature, slave_quadrature,
                  master_facet_number, slave_facet_number, 
                  INTERFACE_SLAVE, 
                  s, num_slaves,
                  L_S);
        assert (!contains_nan(L_S));
        add_global(space, slave_cell_index, slave_dofs, L_S, vec);
      }
    }
  }
}

template < class DataType, int DIM >
template< class LocalAssembler>
void DGGlobalAssembler< DataType, DIM >::assemble_interface_scalar(const VecSpace &space, 
                                                                   LocalAssembler& local_asm, 
                                                                   std::vector< DataType > &values) const 
{

  // Create interface list from mesh
  mesh::ConstMeshPtr mesh = &space.mesh(); // OK since pointer is intrusive

  mesh::InterfaceList if_list = mesh::InterfaceList::create(mesh);

  // Clear and create values data structure
  values.clear();
  values.resize(if_list.size(), 0.);

  DataType L_M, L_S;

  Quadrature< DataType > master_quadrature, slave_quadrature;

  // Loop over interfaces
  size_t i = 0;
  for (mesh::InterfaceList::const_iterator it = if_list.begin(),
       end_it = if_list.end();
       it != end_it; ++it) 
  {
    int remote_index_master = -10;
    mesh->get_attribute_value("_remote_index_", mesh->tdim(),
                              it->master_index(), &remote_index_master);

    // Master dofs
    Element< DataType, DIM > master_elem(space, it->master_index());

    const int master_facet_number = it->master_facet_number();

    L_M = 0.;

    // Initialize master quadrature
    if_q_select_(master_elem, master_elem, 
                 master_facet_number, master_facet_number, 
                 master_quadrature, master_quadrature);

    const int num_slaves = it->num_slaves();
    if (remote_index_master == -1) 
    {
      if (num_slaves == 0) 
      {
        // boundary facet
        local_asm(master_elem, master_elem, 
                  master_quadrature, master_quadrature, 
                  master_facet_number, master_facet_number,
                  INTERFACE_BOUNDARY, 
                  -1, 0,
                  L_M);
      }
      values[i] += L_M;
    }
    // Loop over slaves
    for (int s = 0; s < num_slaves; ++s) 
    {
      int remote_index_slave = -10;
      mesh->get_attribute_value("_remote_index_", mesh->tdim(), it->slave_index(s), &remote_index_slave);
      Element< DataType, DIM > slave_elem(space, it->slave_index(s));
      const int slave_facet_number = it->slave_facet_number(s);

      // Initialize slave quadrature. NB: only called once per slave.
      // default quad selection: all quadrature points lie on that part of the interface,
      // which has non-empty intersection with slave cell
      if_q_select_(master_elem, slave_elem, 
                   master_facet_number, slave_facet_number, 
                   master_quadrature, slave_quadrature);

      if (remote_index_master == -1) 
      {
        // master / slave
        L_S = 0.;
        local_asm(master_elem, slave_elem, 
                  master_quadrature, slave_quadrature,
                  master_facet_number, slave_facet_number, 
                  INTERFACE_SLAVE, 
                  s, num_slaves,
                  L_S);
        values[i] += L_S;
      }
    }
    ++i;
  }
}

template < class DataType, int DIM >
template< class LocalAssembler>
void DGGlobalAssembler< DataType, DIM >::assemble_interface_scalar_cells(const VecSpace &space, 
                                                                         LocalAssembler& local_asm, 
                                                                         std::vector< DataType > &values) const 
{
  // Create interface list from mesh
  mesh::ConstMeshPtr mesh = &space.mesh(); // OK since pointer is intrusive

  mesh::InterfaceList if_list = mesh::InterfaceList::create(mesh);

  // Clear and create values data structure
  values.clear();
  values.resize(mesh->num_entities(mesh->tdim()), 0.);

  int rank = space.dof().my_subdom();

  // Loop over interfaces
  for (mesh::InterfaceList::const_iterator it = if_list.begin(),
       end_it = if_list.end();
       it != end_it; ++it) 
  {
    int remote_index_master = -10;
    mesh->get_attribute_value("_remote_index_", mesh->tdim(),
                              it->master_index(), &remote_index_master);

    // Master dofs
    Element< DataType, DIM > master_elem(space, it->master_index());

    const int master_facet_number = it->master_facet_number();

    DataType L_M, L_S;

    Quadrature< DataType > master_master_quadrature;

    L_M = 0.;

    // Initialize master quadrature
    this->if_q_select_(master_elem, master_elem, master_facet_number,
                       master_facet_number, master_master_quadrature,
                       master_master_quadrature);

    const int num_slaves = it->num_slaves();

    if (remote_index_master == -1) 
    {
      if (num_slaves == 0) 
      {
        // boundary facet
        local_asm(master_elem, master_elem, 
                  master_master_quadrature, master_master_quadrature, 
                  master_facet_number, master_facet_number, 
                  INTERFACE_BOUNDARY, 
                  -1, 0,
                  L_M);
      }

      LOG_DEBUG(3, "[" << rank << "] Master index: " << it->master_index()
                       << " with remote index " << remote_index_master
                       << ", add to master_cell L_MM=" << L_M);

      values[it->master_index()] += L_M;
    }

    // Loop over slaves
    for (int s = 0; s < num_slaves; ++s) 
    {
      int remote_index_slave = -10;
      mesh->get_attribute_value("_remote_index_", mesh->tdim(),
                                it->slave_index(s), &remote_index_slave);
      Element< DataType, DIM > slave_elem(space, it->slave_index(s));
      const int slave_facet_number = it->slave_facet_number(s);

      Quadrature< DataType > master_quadrature, slave_quadrature;

      // Initialize slave quadrature. NB: only called once per slave.
      // default quad selection: all quadrature points lie on that part of the interface,
      // which has non-empty intersection with slave cell
      this->if_q_select_(master_elem, slave_elem, 
                         master_facet_number, slave_facet_number, 
                         master_quadrature, slave_quadrature);

      if (remote_index_master == -1 || remote_index_slave == -1) 
      {
        // master / slave
        L_S = 0.;
        local_asm(master_elem, slave_elem, 
                  master_quadrature, slave_quadrature,
                  master_facet_number, slave_facet_number, 
                  INTERFACE_SLAVE, 
                  s, num_slaves,
                  L_S);

        if (num_slaves > 1 && rank == PRINT_RANK) 
        {
          LOG_DEBUG(2, "[" << rank << "] Master index: " << it->master_index()
                           << " with remote index " << remote_index_master
                           << " and slave index " << it->slave_index(s)
                           << " with remote index " << remote_index_slave
                           << ", add to master_cell L_S=" << 0.5 * L_S);
        }

        // contribution of interface is shared between master and slave cells
        if (remote_index_master == -1) 
        {
          values[it->master_index()] += 0.5 * L_S;
        }
        if (remote_index_slave == -1) 
        {
          values[it->slave_index(s)] += 0.5 * L_S;
        }
      }
    }
  }
}

template < class DataType, int DIM >
template< class LocalAssembler>
void DGGlobalAssembler< DataType, DIM >::assemble_interface_multiple_scalar_cells(const VecSpace &space,
                                                                                  LocalAssembler& local_asm, 
                                                                                  const int num_scalars,
                                                                                  std::vector< typename GlobalAsm::LocalVector > &values) const 
{
  // Create interface list from mesh
  mesh::ConstMeshPtr mesh = &space.mesh(); // OK since pointer is intrusive

  mesh::InterfaceList if_list = mesh::InterfaceList::create(mesh);

  // Clear and create values data structure
  values.clear();
  values.resize(mesh->num_entities(mesh->tdim()));
  for (size_t j = 0; j < values.size(); ++j) 
  {
    values[j].resize(num_scalars, 0.);
  }

#ifndef NDEBUG
  int rank = space.dof().my_subdom();
#endif

  // Loop over interfaces
  for (mesh::InterfaceList::const_iterator it = if_list.begin(), end_it = if_list.end(); it != end_it; ++it) 
  {
    int remote_index_master = -10;
    mesh->get_attribute_value("_remote_index_", mesh->tdim(), it->master_index(), &remote_index_master);

    // Master dofs
    Element< DataType, DIM > master_elem(space, it->master_index());

    const int master_facet_number = it->master_facet_number();

    typename GlobalAsm::LocalVector L_M, L_S;

    Quadrature< DataType > master_master_quadrature;

    // Initialize master quadrature
    this->if_q_select_(master_elem, master_elem, master_facet_number,
                       master_facet_number, master_master_quadrature,
                       master_master_quadrature);

    const int num_slaves = it->num_slaves();

    // Boundary integral
    if (remote_index_master == -1) 
    {
      if (num_slaves == 0) {
        // boundary facet
        local_asm(master_elem, master_elem, 
                  master_master_quadrature, master_master_quadrature, 
                  master_facet_number, master_facet_number, 
                  INTERFACE_BOUNDARY, 
                  -1, 0,
                  L_M);

        LOG_DEBUG(3, "[" << rank << "] Master index: " << it->master_index()
                         << " with remote index " << remote_index_master
                         << ", add to master_cell L_MM="
                         << string_from_range(L_M.begin(), L_M.end()));

        assert(values[it->master_index()].size() == L_M.size());
        for (size_t l = 0; l < L_M.size(); ++l) 
        {
          values[it->master_index()][l] += L_M[l];
        }
      }
    }

    // Interface integrals
    // Loop over slaves
    for (int s = 0; s < num_slaves; ++s) 
    {
      int remote_index_slave = -10;
      mesh->get_attribute_value("_remote_index_", mesh->tdim(), it->slave_index(s), &remote_index_slave);
      Element< DataType, DIM > slave_elem(space, it->slave_index(s));
      const int slave_facet_number = it->slave_facet_number(s);

      Quadrature< DataType > master_quadrature, slave_quadrature;

      // Initialize slave quadrature. NB: only called once per slave.
      // default quad selection: all quadrature points lie on that part of the interface,
      // which has non-empty intersection with slave cell
      this->if_q_select_(master_elem, slave_elem, 
                         master_facet_number, slave_facet_number, 
                         master_quadrature, slave_quadrature);

      if (remote_index_master == -1 || remote_index_slave == -1) 
      {
        local_asm(master_elem, slave_elem, 
                  master_quadrature, slave_quadrature,
                  master_facet_number, slave_facet_number, 
                  INTERFACE_SLAVE, 
                  s, num_slaves, 
                  L_S);

        assert(values[it->master_index()].size() == L_S.size());
        assert(values[it->slave_index(s)].size() == L_S.size());

        // contribution of interface is shared between master and slave cells
        if (remote_index_master == -1) 
        {
          for (size_t l = 0; l < L_S.size(); ++l) 
          {
            values[it->master_index()][l] += 0.5 * L_S[l];
          }
        }
        if (remote_index_slave == -1) 
        {
          for (size_t l = 0; l < L_S.size(); ++l) 
          {
            values[it->slave_index(s)][l] += 0.5 * L_S[l];
          }
        }
      }
    }
  }
}

template < class DataType, int DIM >
void DGGlobalAssembler< DataType, DIM >::distribute_interface_to_cell_values_naive(const VecSpace &space, 
                                                                                   std::vector< DataType > &cell_values,
                                                                                   const std::vector< DataType > &interface_values) const 
{
  // Create interface list from mesh
  mesh::ConstMeshPtr mesh = &space.mesh(); // OK since pointer is intrusive

  // Check compatibility of mesh and cell_values vector
  assert(mesh->num_entities(mesh->tdim()) == cell_values.size());

  mesh::InterfaceList if_list = mesh::InterfaceList::create(mesh);

  // Check compatibility of interface list and interface_values vector
  assert(if_list.size() == interface_values.size());

  // Loop over interfaces
  int i = 0;
  for (mesh::InterfaceList::const_iterator it = if_list.begin(),
                                     end_it = if_list.end();
       it != end_it; ++it) {
    int remote_index_master = -10;
    mesh->get_attribute_value("_remote_index_", mesh->tdim(),
                              it->master_index(), &remote_index_master);

    const int num_slaves = it->num_slaves();
    if (remote_index_master == -1) {
      if (num_slaves > 0) {
        // Master only gets half of the contribution of interface value
        cell_values[it->master_index()] += 0.5 * interface_values[i];
      } else {
        // boundary facet
        // Master only gets contribution of interface value
        cell_values[it->master_index()] += interface_values[i];
      }
    }

    // weight per slave
    DataType weight_slave = 0.5;
    if (num_slaves > 0) {
      weight_slave /= num_slaves;
    }
    // Loop over slaves
    for (int s = 0; s < num_slaves; ++s) {
      int remote_index_slave = -10;
      mesh->get_attribute_value("_remote_index_", mesh->tdim(),
                                it->slave_index(s), &remote_index_slave);

      if (remote_index_slave == -1) {
        cell_values[it->slave_index(s)] += weight_slave * interface_values[i];
      }
    }
    ++i;
  }
}

template < class DataType, int DIM >
void DGGlobalAssembler< DataType, DIM >::set_interface_quadrature_selection_fun( IFQuadratureSelectionFun q_select) 
{
  if_q_select_ = q_select;
}

} // namespace hiflow

#endif
