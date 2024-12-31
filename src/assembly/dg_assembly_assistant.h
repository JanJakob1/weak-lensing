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

#ifndef _DG_ASSEMBLY_ASSISTANT_H_
#define _DG_ASSEMBLY_ASSISTANT_H_

#include "assembly/assembly_assistant.h"
#include "assembly/dg_assembly.h"
#include "assert.h"

namespace hiflow {
using doffem::FEType;

/// \author Staffan Ronnas, Jonathan Schwegler, Simon Gawlok, Philipp Gerstner

/// \brief Provides assembly functionalities needed for DG, i.e., local
/// assembly over interfaces where function values are needed from
/// several adjacent cells.

template < int DIM, class DataType > class DGAssemblyAssistant {
public:
  typedef typename DGGlobalAssembler< DataType, DIM >::InterfaceSide InterfaceSide;
  typedef la::SeqDenseMatrix< DataType > LocalMatrix;
  typedef std::vector< DataType > LocalVector;

  /// \brief Standard constructor
  DGAssemblyAssistant();

  /// \brief Get AssemblyAssistant corresponding to trial space finite element

  inline AssemblyAssistant< DIM, DataType > &trial() {
    assert(this->trial_aa_ != nullptr);
    return *(this->trial_aa_);
  }
  
  inline const AssemblyAssistant< DIM, DataType > &trial() const {
    assert(this->trial_aa_ != nullptr);
    return *(this->trial_aa_);
  }
  
  /// \brief Get AssemblyAssistant corresponding to test space finite element
  inline AssemblyAssistant< DIM, DataType > &test() {
    assert(this->test_aa_ != nullptr);
    return *(this->test_aa_);
  }

  inline const AssemblyAssistant< DIM, DataType > &test() const {
    assert(this->test_aa_ != nullptr);
    return *(this->test_aa_);
  }
  
  inline AssemblyAssistant< DIM, DataType >& master() {
    return (this->master_aa_);
  }

  inline const AssemblyAssistant< DIM, DataType >& master() const {
    return (this->master_aa_);
  }  
  
  inline AssemblyAssistant< DIM, DataType >& slave() {
    // slave AA is not set up on boundary interfaces
    assert (!this->is_boundary_);
    return (this->slave_aa_);
  }

  inline const AssemblyAssistant< DIM, DataType >& slave() const {
    // slave AA is not set up on boundary interfaces
    assert (!this->is_boundary_);
    return (this->slave_aa_);
  }
     
  /// \brief Get number of variables
  /// Provided by test space finite element

  inline int num_vars() const {
    return this->master_aa_.num_vars();
  }

  inline int num_fe() const {
    return this->master_aa_.num_fe();
  }
  
  // AssemblyAssistant functions that are the same for trial and test function
  /// \brief Get number of quadrature points on current interface;
  /// Provided by test space finite element
  inline int num_quadrature_points() const {
    return this->master_aa_.num_quadrature_points();
  }

  inline int num_dofs_total() const {
    return this->master_aa_.num_dofs_total();
  }
    
  /// \brief Get coordinates of a quadrature point on current interface;
  /// Provided by test space finite element
  /// \param q Number of quadrature point
  inline Vec< DIM, DataType > x(int q) const {
    return this->master_aa_.x(q);
  }
  
  /// \brief Get quadrature weight at a quadrature point on current interface;
  /// Provided by test space finite element
  /// \param q Number of quadrature point
  inline DataType w(int q) const {
    return this->master_aa_.w(q);
  }
  
  /// \brief Get surface element at a quadrature point on current interface;
  /// Provided by test space finite element
  /// \param q Number of quadrature point
  inline DataType ds(int q) const {
    return this->ds_aa_->ds(q);
  }
  
  inline DataType detJ(int q) const {
    return this->master_aa_.detJ(q);
  }
  
  inline size_t first_dof_for_var (size_t var) const {
    return this->master_aa_.first_dof_for_var(var);
  }

  inline size_t last_dof_for_var (size_t var) const {
    return this->master_aa_.last_dof_for_var(var);
  }

  inline size_t num_dofs (size_t var) const 
  {
    return this->master().num_dofs_for_fe(this->master().var_2_fe (var));
  }

  inline size_t dof_index(size_t i, size_t var) const 
  {
    return this->master().dof_offset_for_var(var) + i;
  }

  inline DataType phi(int s, int q, int var) const 
  {
    return this->master().phi(s,q,var);
  }

  inline DataType Phi(int i, int q, int var) const 
  {
    return this->master().Phi(i,q,var);
  }

  inline const Vec< DIM, DataType > & grad_phi(int s, int q, int var) const 
  {
    return this->master().grad_phi(s,q,var);
  }

  inline const Vec< DIM, DataType > & grad_Phi(int i, int q, int var) const 
  {
    return this->master().grad_Phi(i,q,var);
  }

  inline const Mat< DIM, DIM, DataType > & H_phi(int s, int q, int var) const 
  {
    return this->master().H_phi(s, q, var);
  }

  inline const Mat< DIM, DIM, DataType > & H_Phi(int s, int q, int var) const 
  {
    return this->master().H_Phi(s, q, var);
  }

  template < class VectorType >
  void evaluate_fe_function(const VectorType &coefficients, 
                            int var,
                            FunctionValues< DataType > &function_values) const
  {
    this->master().evaluate_fe_function(coefficients, var, function_values);
  }

  template < class VectorType >
  void evaluate_fe_function_gradients(const VectorType &coefficients, 
                                      int var,
                                      FunctionValues< Vec< DIM, DataType > > &function_gradients) const
  {
    this->master().evaluate_fe_function_gradients(coefficients, var, function_gradients);
  }

  template < class VectorType >
  void evaluate_fe_function_hessians(const VectorType &coefficients, 
                                     int var,
                                     FunctionValues< Mat< DIM, DIM, DataType > > &function_hessians) const
  {
    this->master().evaluate_fe_function_hessians(coefficients, var, function_hessians);
  }
                                     
  void initialize_for_element(const Element< DataType, DIM > &element,
                              const Quadrature< DataType > &element_quadrature,
                              bool compute_hessians);
 
  void initialize_for_facet(const Element< DataType, DIM > &element,
                            const Quadrature< DataType > &facet_quadrature,
                            int facet_number, 
                            bool compute_hessians);

  /// \brief Initialize an interface between two elements (trial_elem,
  /// test_elem). The physical points of the two quadratures must be the same.
  /// \param master_elem master space finite element
  /// \param slave_elem slave space finite element
  /// \param master_quad Quadrature on master space finite element
  /// \param slave_quad Quadrature on slave space finite element
  /// \param master_facet_number Facet number of master space finite element
  /// \param slave_facet_number Facet number of slave space finite element
  /// \param trial_if_side Kind of interface of left space finite element
  /// \param test_if_side Kind of interface of right space finite element
  void initialize_for_interface(const Element< DataType, DIM > &master_elem,
                                const Element< DataType, DIM > &slave_elem,
                                const Quadrature< DataType > &master_quad,
                                const Quadrature< DataType > &slave_quad,
                                int master_facet_number, 
                                int slave_facet_number,
                                InterfaceSide trial_if_side,
                                InterfaceSide test_if_side,
                                bool compute_hessians);
                                
protected:
  /// AssemblyAssistant corresponding to "master" finite element
  AssemblyAssistant< DIM, DataType > master_aa_;
  /// AssemblyAssistant corresponding to "slave" finite element
  AssemblyAssistant< DIM, DataType > slave_aa_;
  /// AssemblyAssistant corresponding to trial space finite element
  AssemblyAssistant< DIM, DataType > *trial_aa_;
  /// AssemblyAssistant corresponding to test space finite element
  AssemblyAssistant< DIM, DataType > *test_aa_;

  AssemblyAssistant< DIM, DataType > *ds_aa_;

  /// Flag if current interface is a boundary facet
  bool is_boundary_;

  /// Interface sides of neighbouring cells
  InterfaceSide trial_if_side_, test_if_side_;
};

template < int DIM, class DataType >
DGAssemblyAssistant< DIM, DataType >::DGAssemblyAssistant() {}

template < int DIM, class DataType >
void DGAssemblyAssistant< DIM, DataType >::initialize_for_element(const Element< DataType, DIM > &element,
                                                                  const Quadrature< DataType > &element_quadrature,
                                                                  bool compute_hessians)
{
  this->master_aa_.initialize_for_element(element, element_quadrature, compute_hessians);
  this->test_aa_     = &(this->master_aa_);
  this->trial_aa_    = &(this->master_aa_);
  this->ds_aa_       = &(this->master_aa_);
  this->is_boundary_ = false;
  this->trial_if_side_ = DGGlobalAssembler< DataType, DIM >::INTERFACE_NONE;
  this->test_if_side_  = DGGlobalAssembler< DataType, DIM >::INTERFACE_NONE;
}
                              
template < int DIM, class DataType >
void DGAssemblyAssistant< DIM, DataType >::initialize_for_facet(const Element< DataType, DIM > &element,
                                                                const Quadrature< DataType > &facet_quadrature,
                                                                int facet_number, 
                                                                bool compute_hessians)
{
  this->master_aa_.initialize_for_facet(element,
                                       facet_quadrature,
                                       facet_number, 
                                       compute_hessians);
                                      
  this->test_aa_     = &(this->master_aa_);
  this->trial_aa_    = &(this->master_aa_);
  this->ds_aa_       = &(this->master_aa_);
  // TODO: actually, we don't know this, right?
  this->is_boundary_ = false;
  this->trial_if_side_ = DGGlobalAssembler< DataType, DIM >::INTERFACE_NONE;
  this->test_if_side_  = DGGlobalAssembler< DataType, DIM >::INTERFACE_NONE;
}
                            
template < int DIM, class DataType >
void DGAssemblyAssistant< DIM, DataType >::initialize_for_interface(const Element< DataType, DIM > &master_elem, 
                                                                    const Element< DataType, DIM > &slave_elem,
                                                                    const Quadrature< DataType > &master_quad,
                                                                    const Quadrature< DataType > &slave_quad, 
                                                                    int master_facet_number,
                                                                    int slave_facet_number, 
                                                                    InterfaceSide trial_if_side,
                                                                    InterfaceSide test_if_side,
                                                                    bool compute_hessians) 
{
  // Convention: left -- trial  functions, right -- test functions. 
  this->trial_if_side_ = trial_if_side;
  this->test_if_side_  = test_if_side;
  this->is_boundary_   = (test_if_side == DGGlobalAssembler< DataType, DIM >::INTERFACE_BOUNDARY);
  
  // initialize master / slave assembly assistants
  this->master_aa_.initialize_for_facet(master_elem, master_quad, master_facet_number, compute_hessians);
  
  if (!(this->is_boundary_)) 
  {
    this->slave_aa_.initialize_for_facet(slave_elem, slave_quad, slave_facet_number, compute_hessians);  
    this->ds_aa_ = &(this->slave_aa_);
  }
  else
  {
    this->ds_aa_ = &(this->master_aa_);
  }
    
  // set trial / test assembly assistant pointers to master / slave
  if (this->is_boundary_) 
  { 
    this->test_aa_  = &(this->master_aa_);
    this->trial_aa_ = &(this->master_aa_);
  } 
  else 
  {
    switch (trial_if_side_)
    {
      case DGGlobalAssembler< DataType, DIM >::INTERFACE_MASTER:
        this->trial_aa_ = &(this->master_aa_);
        break;
      case DGGlobalAssembler< DataType, DIM >::INTERFACE_SLAVE:
        this->trial_aa_ = &(this->slave_aa_);
        break;
      case DGGlobalAssembler< DataType, DIM >::INTERFACE_NONE:
        this->trial_aa_ = &(this->master_aa_);
        break;
      default:
        assert (false);
        break;
    } 
    switch (test_if_side_)
    {
      case DGGlobalAssembler< DataType, DIM >::INTERFACE_MASTER:
        this->test_aa_ = &(this->master_aa_);
        break;
      case DGGlobalAssembler< DataType, DIM >::INTERFACE_SLAVE:
        this->test_aa_ = &(this->slave_aa_);
        break;
      default:
        assert(false);
        break;
    }
  }
  assert(this->trial_aa_->num_quadrature_points() == this->test_aa_->num_quadrature_points());
}


/*
 * OLD IMPLEMENTATION
{
  this->trial_if_side_ = trial_if_side;
  this->test_if_side_ = test_if_side;
  this->is_boundary_ =
      (test_if_side == DGGlobalAssembler< DataType, DIM >::INTERFACE_BOUNDARY);
      

  this->master_aa_.initialize_for_facet(test_elem, test_quad,
                                        test_facet_number, compute_hessians);
  this->test_aa_ = &(this->master_aa_);
  
  // only init both facets if they are different.
  
  if (this->is_boundary_ || trial_if_side == test_if_side) 
  { // if trial elem == test elem
    this->trial_aa_ = this->test_aa_;
  } 
  else 
  {
    this->slave_aa_.initialize_for_facet(trial_elem, trial_quad, trial_facet_number, compute_hessians);
    this->trial_aa_ = &(this->slave_aa_);
  }

  assert(this->trial_aa_->num_quadrature_points() == this->test_aa_->num_quadrature_points());
}
*/
} // namespace hiflow
#endif
