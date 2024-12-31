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

#ifndef _ASSEMBLY_ASSISTANT_H_
#define _ASSEMBLY_ASSISTANT_H_

#include <algorithm>
#include <numeric>
#include <cassert>
#include <vector>

#include "assembly/assembly_assistant_values.h"
#include "assembly/function_values.h"
#include "common/pointers.h"
#include "common/vector_algebra.h"
#include "dof/dof_fem_types.h"
#include "fem/cell_trafo/cell_transformation.h"
#include "fem/fe_reference.h"
#include "fem/reference_cell.h"
#include "linear_algebra/seq_dense_matrix.h"
#include "mesh/cell_type.h"
#include "mesh/types.h"
#include "quadrature/quadrature.h"
#include "space/element.h"

///
/// \file assembly_assistant.h Convenience functionality for local assembly.
/// \author Staffan Ronnas<br>Simon Gawlok<br> Philipp Gerstner
/// \see concept_assembly

namespace hiflow {

// TODO: if this type of assembly is desired, talk to Michael about
// changing Quadrature, etc to use Vec & Mat classes

/// \brief Helper class for local assembly and integration functors.
/// \see concept_assembly

template < int DIM, class DataType > 
class AssemblyAssistant {
public:
  // Typenames for local assembly
  typedef la::SeqDenseMatrix< DataType > LocalMatrix;
  typedef std::vector< DataType > LocalVector;

  // TODO: let user decide what to compute through parameter
  AssemblyAssistant();

  // TODO: move these functions outside class and put documentation in here.

  /// \brief Shape function values on physical element.
  ///
  /// \param s    index of the shape function (relative to variable var if given) 
  /// \param q    index of the quadrature point 
  /// \param var  variable number 
  /// \return     Shape function value \f$\varphi_{v,i}(\xi_q)\f$
  inline DataType phi(int s, int q, int var = 0) const;

  /// \brief Shape function values on physical element.
  ///
  /// \param i    index of the shape function (absolute on given element) 
  /// \param q    index of the quadrature point 
  /// \param var  variable number 
  /// \return     Shape function value \f$\varphi_{v,i}(\xi_q)\f$
  inline DataType Phi(int i, int q, int var) const;
  
  /// \brief Shape function gradients on physical element.
  ///
  /// \param s    index of the shape function (relative to variable var if given) 
  /// \param q    index of the quadrature point 
  /// \param var  variable number 
  /// \return     Shape function gradient \f$\nabla\varphi_{v,i}(x_q)\f$
  inline const Vec< DIM, DataType > &grad_phi(int s, int q, int var = 0) const;
  
  inline const Vec< DIM, DataType > &grad_Phi(int i, int q, int var) const;

  /// \brief Shape function hessian on physical element.
  ///
  /// \param s    index of the shape function (relative to variable var if given) 
  /// \param q    index of the quadrature point 
  /// \param var  variable number 
  /// \return     Shape function hessian
  /// \f$\partial_j\partial_k\varphi_{v,i}(x_q)\f$
  ///     odo   This function has not been tested yet.
  inline const Mat< DIM, DIM, DataType > &H_phi(int s, int q, int var = 0) const;

  inline const Mat< DIM, DIM, DataType > &H_Phi(int i, int q, int var) const;

  /// \brief Shape function laplacian on physical element.
  ///
  /// \param s    index of the shape function (relative to variable var if given) 
  /// \param q    index of the quadrature point 
  /// \param var  variable number 
  /// \return     Shape function laplacian
  /// \f$\partial_j\partial_k\varphi_{v,i}(x_q)\f$
  ///     odo   This function has not been tested yet.
  inline DataType L_phi(int s, int q, int var = 0) const;

  inline DataType L_Phi(int i, int q, int var) const;

  /// \brief Quadrature point on physical element.
  /// \param q    index of the quadrature point
  /// \return     Quadrature point \f$x_q\f$ on physical element
  inline const Vec< DIM, DataType > &x(int q) const;

  /// \brief Jacobian matrix of cell transformation.
  /// \param q    index of the quadrature point
  /// \return     Jacobian matrix \f$DF_K(\xi_q)\f$ of the cell transformation
  /// \f$F_K\f$.
  inline const Mat< DIM, DIM, DataType > &J(int q) const;

  /// \brief Inverse transpose of Jacobian matrix of cell transformation.
  /// \param q    index of the quadrature point
  /// \return     Inverse transpose \f$(DF_K(\xi_q))^{-T}\f$ of the jacobian
  /// matrix of the cell transformation \f$F_K\f$.
  inline const Mat< DIM, DIM, DataType > &JinvT(int q) const;

  /// \brief Hessian tensor of cell transformation.
  ///
  /// \param q    index of the quadrature point
  ///
  /// \return Hessian tensor \f$H_{ijk} =
  /// \partial_j\partial_k{F_i}\f$ of the cell transformation
  /// \f$F_K\f$. Index i is the position in the vector, and jk
  /// the indices for the contained matrices.
  inline const std::vector< Mat< DIM, DIM, DataType > > &H_F(int q) const;

  /// \brief Determinant of Jacobian matrix of cell transformation.
  ///
  /// \param q    index of the quadrature point
  /// \return     Jacobian determinant \f$det(DF_K(\xi_q))\f$ of the cell
  /// transformation \f$F_K\f$.
  inline DataType detJ(int q) const;

  /// \brief Jacobian matrix projected onto current facet.
  ///
  /// \pre        Assembly assistant was initialized on a facet.
  /// \param q    index of the quadrature point
  /// \return     Jf(q) = J(q) * R, where R is projection matrix corresponding
  /// to current facet.
  inline const Mat< DIM, DIM - 1, DataType > &Jf(int q) const;

  /// \brief Surface integration element.
  ///
  /// \pre        Assembly assistant was initialized on a facet.
  /// \param q    index of the quadrature point
  /// \return     Surface integration element ds(q) = \sqrt(det(Jf^T * Jf)) on
  /// current facet.
  inline DataType ds(int q) const;

  /// \brief Surface normal.
  ///
  /// \pre        Assembly assistant was initialized on a facet.
  /// \param q    index of the quadrature point
  /// \return     Surface normal n on current (reference) facet.
  inline const Vec< DIM, DataType > &n(int q) const;

  /// \brief The Nitsche regularization parameter used in
  /// Nitsche method.
  ///
  /// \pre        Assembly assistant was initialized on a facet.
  /// \param q    index of the quadrature point
  /// \return     Nitsche regularization parameter.
  ///     odo       Compute it in some good way.
  inline DataType nitsche_regularization(int q) const;

  /// \brief Quadrature point on reference element.
  ///
  /// \param q    index of the quadrature point
  /// \return     Quadrature point \f$\xi_q\f$ on reference element
  inline const Vec< DIM, DataType > &q_point(int q) const;

  /// \brief Number of quadrature points.
  ///
  /// \return Number of quadrature points for current quadrature rule.
  inline int num_quadrature_points() const;

  /// \brief Quadrature weights.
  ///
  /// \param q    index of the quadrature point
  /// \return     Quadrature weight \f$w_q\f$ for the
  inline DataType w(int q) const;

  /// \brief Number of different FE ansatz spaces.
  /// Note: 0 <= fe index < nb_fe
  inline size_t nb_fe() const;

  /// \brief number of physical variables, >= nb_fe
  inline size_t nb_var() const;

  /// \brief number of different fe types. 
  /// Note: nb_fe_types >= nb_fe with nb_fe_types > nb_fe in case of p-refinement 
  /// and/or mesh consisting of different types of cells
  inline size_t nb_fe_types() const; 
  
  /// \brief return index in RefElement array for given fe index 
  /// (relevant in case of varying cell types and/or polynomial degrees within a single fe ansatz) 
  inline size_t fe_2_fe_type (size_t fe) const; 

  /// \brief return index in RefElement array for given physical variable
  /// (relevant in case of varying cell types and/or polynomial degrees within a single fe ansatz) 
  inline size_t var_2_fe_type (size_t var) const; 

  /// \brief return fe index for given physical variable
  inline size_t var_2_fe (size_t var) const; 

  /// \brief return fe component of given physical variable
  /// Relevant in case of vector valued FEs. For scalar valued FEs: comp = 0
  inline size_t var_2_comp (size_t var) const; 

  /// \brief return first basis function index with nonzero entry of given physical variable
  inline size_t first_dof_for_var (size_t var) const; 

  /// \brief return last+1 basis function index with nonzero entry of given physical variable
  inline size_t last_dof_for_var (size_t var) const; 
    
  inline size_t dof_offset_for_fe (size_t fe) const; 

  inline size_t dof_offset_for_var (size_t var) const; 

  /// \return Number of local dofs for variable var.
  inline size_t num_dofs (size_t var) const; 

  inline size_t num_dofs_for_fe (size_t fe) const; 

  /// \brief Number of local dofs.
  inline size_t num_dofs_total() const; 

  /// \brief Local dof index for dof s of variable var.
  /// \param s     dof index relative to variable var ( in range [0, num_dofs(var)) ) 
  /// \param var   variable number \return Local dof index \f$l = l(var, s)\f$ for dof \f$s\f$ for variable var.
  inline size_t dof_index(size_t i, size_t var) const;

  inline size_t dof_index_for_fe(size_t i, size_t fe) const; 

  /// \brief Element size parameter h
  /// \return Element size parameter h
  inline DataType h() const;

  bool need_cell_trafo_hessian (bool compute_hessians) const;
  bool need_cell_trafo_grad_inv_detJ (bool compute_hessians) const;

  /// \brief Evaluate a finite element function on the current element.
  ///
  /// \param[in] coefficients  global vector of dof values for the function.
  /// \param[in] var           variable to compute the function for.
  /// \param[out] function_values  function values evaluated at all quadrature points.
  template < class VectorType >
  void evaluate_fe_function(const VectorType &coefficients, 
                            int var,
                            FunctionValues< DataType > &function_values) const;

  /// \brief Evaluate the gradient of a finite element function on the current element.

  /// \param[in] coefficients  global vector of dof values for the function.
  /// \param[in] var           variable to compute the gradient for.
  /// \param[out] function_gradients  function gradients at all quadrature points.
  template < class VectorType >
  void evaluate_fe_function_gradients(const VectorType &coefficients, 
                                      int var,
                                      FunctionValues< Vec< DIM, DataType > > &function_gradients) const;

  /// \brief Evaluate the hessian of a finite element function on the current element.
  ///
  /// \param[in] coefficients  global vector of dof values for the function.
  /// \param[in] var           variable to compute the gradient for.
  /// \param[out] function_hessians  function hessians at all quadrature points.
  template < class VectorType >
  void evaluate_fe_function_hessians(const VectorType &coefficients, 
                                     int var,
                                     FunctionValues< Mat< DIM, DIM, DataType > > &function_hessians) const;

  /// \brief Evaluate the hessian and laplacian of a finite element function on the current element.
  ///
  /// \param[in] coefficients  global vector of dof values for the function.
  /// \param[in] var           variable to compute the gradient for.
  /// \param[out] function_hessians  function hessians at all quadrature points.
  /// \param[out] function_laplace  function laplacians at all quadrature points.

  template < class VectorType >
  void evaluate_fe_function_hessians_and_laplace( const VectorType &coefficients, 
                                                  int var,
                                                  FunctionValues< Mat< DIM, DIM, DataType > > &function_hessians,
                                                  FunctionValues< DataType > &function_laplace) const;

  /// \brief Initialize the assistant for an element and a quadrature.
  ///
  /// \details Recomputes all necessary values on the given element for the given quadrature formula
  /// \param[in] element     element for which to initialize the AssemblyAssistant 
  /// \param[in] quadrature  quadrature rule
  void initialize_for_element(const Element< DataType, DIM > &element,
                              const Quadrature< DataType > &element_quadrature,
                              bool compute_hessians);

  /// \brief Initialize the assistant for a facet of the element and a
  /// quadrature.
  ///
  /// \details Recomputes all necessary values on the given
  /// element for the given quadrature formula and facet.
  ///
  /// \param[in] element element for which to initialize the AssemblyAssistant 
  /// \param[in] facet_quadrature quadrature rule (should have points only on the given element facet)
  /// \param[in] facet_number local number of facet in element
  
  void initialize_for_facet(const Element< DataType, DIM > &element,
                            const Quadrature< DataType > &facet_quadrature,
                            int facet_number, 
                            bool compute_hessians);

private:
  // initialization sub-functions
  void initialize( const Element< DataType, DIM > &element,
                   const Quadrature< DataType > &quadrature,
                   int facet_number,
                   bool is_facet,
                   bool compute_hessians);
                                                     
  bool update_ref_quadrature(const Quadrature< DataType > &new_quadrature,
                             bool force_update);
                             
  bool update_fe_types(const Element< DataType, DIM > &element,
                       bool force_update);
  
  bool update_cell_trafo(const Element< DataType, DIM > &element,
                         bool changed_quad,
                         bool compute_hessians);

  void compute_ref_cell_values(bool compute_hessians);
  void compute_trans_cell_values(bool compute_hessians);
  void compute_trans_facet_values();
  void compute_facet_projection_matrix();
  void compute_facet_normal();

  int find_fe_type(doffem::RefElement< DataType, DIM > const * fe_type) const;
  
  template < class VectorType >
  void extract_dof_values(const VectorType &coefficients, 
                          size_t fe_type,
                          std::vector<DataType>& local_coeff) const;

  /// Current cell index
  mesh::EntityNumber cell_index_; 
  
  /// Current facet number (-1 on cell)
  int facet_number_;              
  
  /// Current facet projection matrix
  Mat< DIM, DIM - 1, DataType > facet_proj_;        
  
  /// Current facet normal
  Vec< DIM, DataType > n_; 
  
  /// Normals on physical element
  FunctionValues< Vec< DIM, DataType > > mapped_n_; 
  
  /// Regularization parameter for Nitsche method
  DataType nitsche_regularization_; 
  
  /// Current quadrature
  Quadrature< DataType > quadrature_;                  
  
  /// Current FEType:s
  std::vector< doffem::RefElement< DataType, DIM > const * > fe_types_; 
  
  /// Pointer to check if space changes
  const VectorSpace< DataType, DIM > *space_;  

  /// number of physical variables
  size_t num_var_;
  
  /// Mapping var -> comp in RefElement
  std::vector< size_t > var_2_comp_; 
  
  /// Mapping RefElement in complete FE VectorSpace -> unique RefElement in list fe_types_
  std::vector< size_t > fe_2_fe_type_; 

  /// Mapping var -> RefElement in complete FE VectorSpace
  std::vector< size_t > var_2_fe_; 
  
  /// Mapping var ->unique RefElement in list fe_types_
  std::vector< size_t > var_2_fe_type_;

  /// mapping var -> fist and last basis index with nonzero component var
  std::vector< size_t > first_dof_for_var_;
  std::vector< size_t > last_dof_for_var_;
  
  /// Number of dofs per RefElement in complete FE VectorSpace
  std::vector< size_t > num_dofs_for_fe_;   
  
  /// total number of dofs on current element
  size_t num_dofs_;

  /// Dof offsets per RefElement in complete FE VectorSpace
  std::vector< size_t > dof_offsets_for_fe_;
  std::vector< size_t > phi_offsets_for_var_;
   
  /// Global dofs for current element
  std::vector< int > global_dof_indices_; 

  std::vector< DataType > dof_factors_;
  
  /// Cell transform for current element
  std::shared_ptr<doffem::CellTransformation< DataType, DIM > const> cell_transform_;
  
  /// Cell type of current element           
  doffem::ConstRefCellPtr<DataType, DIM> ref_cell_;
  
  /// Element size parameter h
  DataType h_;                    

  /// Quadrature points on reference element
  std::vector< Vec< DIM, DataType > > q_points_;          
  std::vector< Vec< DIM, DataType > > q_points_prev_;    
  
  /// Quadrature weights          
  std::vector< DataType > weights_; 

  /// indexing of the following variables:
  /// psi[q][fe][comp * dim_fe + j] :
  /// q: index of quad point
  /// fe: index of RefElement
  /// dim_fe: number of basis functions within RefElement fe 
  /// comp: index of component within RefElement fe
  /// j: index of basis function within RefElement fe
  
  /// Shape function values on reference element
  FunctionValues< std::vector< std::vector< DataType > > > phi_hat_; 

  /// Shape function values
  FunctionValues< std::vector< std::vector< DataType > > > fe_phi_; 
  FunctionValues< std::vector< DataType > > phi_;

  /// Shape function gradients on reference element
  FunctionValues< std::vector< std::vector< Vec< DIM, DataType > > > > grad_phi_hat_;

  /// Shape function hessians on reference element
  FunctionValues< std::vector< std::vector< Mat< DIM, DIM, DataType > > > > H_phi_hat_;

  /// Shape function gradients on physical element
  FunctionValues< std::vector<std::vector< Vec< DIM, DataType > > > > fe_grad_phi_;
  FunctionValues< std::vector< Vec< DIM, DataType > > > grad_phi_;
  
  /// Shape function hessians on physical element
  FunctionValues< std::vector<std::vector< Mat< DIM, DIM, DataType > > > > fe_H_phi_;
  FunctionValues< std::vector< Mat< DIM, DIM, DataType > > > H_phi_;
  
  /// Shape function Laplacian 
  FunctionValues< std::vector< DataType > > L_phi_;

  /// H-mapped shape function gradients on physical element
  FunctionValues< std::vector<std::vector< Mat< DIM, DIM, DataType > > > > H_mapped_grad_;
  
  /// Quadrature points on physical element
  FunctionValues< Vec< DIM, DataType > > x_;
   
  /// Jacobian matrix of cell transformation at quadrature points
  FunctionValues< Mat< DIM, DIM, DataType > > J_;
  
  /// Determinants of jacobian matrices
  FunctionValues< DataType > detJ_; 
  
  /// Gradient of inverse determinants of jacobian matrices
  FunctionValues< Vec<DIM,DataType > > grad_inv_detJ_; 
  
  /// Inverse-transpose of jacobian matrices
  FunctionValues< Mat< DIM, DIM, DataType > > Jinv_; 
  
  /// Inverse-transpose of jacobian matrices
  FunctionValues< Mat< DIM, DIM, DataType > > JinvT_; 
  
  /// Jacobian matrix projected onto facet
  FunctionValues< Mat< DIM, DIM - 1, DataType > > Jf_;

  ///< Hessian of cell transformation at quadrature points
  FunctionValues< std::vector< Mat< DIM, DIM, DataType > > > H_;
    
  /// Surface element on facet
  FunctionValues< DataType > ds_; 
  
  Vec<DIM, DataType> zero_vec_;
  Mat<DIM, DIM, DataType> zero_mat_;

  
  bool initialized_fe_types_;
  
  mutable std::vector< DataType > local_coefficients_;
  mutable std::vector< int > tmp_indices_;
  mutable std::vector< DataType > tmp_factors_;
  
};

template < int DIM, class DataType >
AssemblyAssistant< DIM, DataType >::AssemblyAssistant() 
: cell_transform_(0), space_(nullptr)
{
  assert(DIM > 0);
  assert(DIM <= 3);
  facet_number_ = -1;
  num_var_ = -1;
  cell_index_ = -1;
  this->initialized_fe_types_ = false;
}

template < int DIM, class DataType >
inline size_t AssemblyAssistant< DIM, DataType >::fe_2_fe_type (size_t fe) const 
{
  assert (fe < fe_2_fe_type_.size());
  return this->fe_2_fe_type_[fe];
}

template < int DIM, class DataType >
inline size_t AssemblyAssistant< DIM, DataType >::var_2_fe_type (size_t var) const 
{
  assert (var < var_2_fe_type_.size());
  return this->var_2_fe_type_[var];
}

template < int DIM, class DataType >
inline size_t AssemblyAssistant< DIM, DataType >::var_2_fe (size_t var) const 
{
  assert (var < var_2_fe_.size());
  return this->var_2_fe_[var];
}

template < int DIM, class DataType >
inline size_t AssemblyAssistant< DIM, DataType >::var_2_comp (size_t var) const 
{
  assert (var < var_2_comp_.size());
  return this->var_2_comp_[var];
}

template < int DIM, class DataType >
inline size_t AssemblyAssistant< DIM, DataType >::first_dof_for_var (size_t var) const 
{
  assert (var < first_dof_for_var_.size());
  return this->first_dof_for_var_[var];
}

template < int DIM, class DataType >
inline size_t AssemblyAssistant< DIM, DataType >::last_dof_for_var (size_t var) const 
{
  assert (var < last_dof_for_var_.size());
  return this->last_dof_for_var_[var];
}

template < int DIM, class DataType >
inline size_t AssemblyAssistant< DIM, DataType >::dof_offset_for_fe (size_t fe) const 
{
  assert (fe < dof_offsets_for_fe_.size());
  return this->dof_offsets_for_fe_[fe];
}

template < int DIM, class DataType >
inline size_t AssemblyAssistant< DIM, DataType >::dof_offset_for_var (size_t var) const 
{
  return this->dof_offset_for_fe( this->var_2_fe (var) );
}

template < int DIM, class DataType >
inline size_t AssemblyAssistant< DIM, DataType >::num_dofs_total () const 
{
  return this->num_dofs_;
}

template < int DIM, class DataType >
inline size_t AssemblyAssistant< DIM, DataType >::num_dofs (size_t var) const 
{
  return this->num_dofs_for_fe(this->var_2_fe (var));
}

template < int DIM, class DataType >
inline size_t AssemblyAssistant< DIM, DataType >::num_dofs_for_fe (size_t fe) const 
{
  assert (fe < num_dofs_for_fe_.size());
  return this->num_dofs_for_fe_[fe];
}

template < int DIM, class DataType >
inline size_t AssemblyAssistant< DIM, DataType >::nb_var() const 
{
  return num_var_;
}

template < int DIM, class DataType >
inline size_t AssemblyAssistant< DIM, DataType >::nb_fe() const 
{
  return num_dofs_for_fe_.size();
}

template < int DIM, class DataType >
inline size_t AssemblyAssistant< DIM, DataType >::nb_fe_types() const 
{
  return fe_types_.size();
}

template < int DIM, class DataType >
inline size_t AssemblyAssistant< DIM, DataType >::dof_index(size_t i, size_t var) const 
{
  return this->dof_offset_for_var(var) + i;
}

template < int DIM, class DataType >
inline size_t AssemblyAssistant< DIM, DataType >::dof_index_for_fe(size_t i, size_t fe) const 
{
  return this->dof_offset_for_fe(fe) + i;
}

template < int DIM, class DataType >
inline DataType AssemblyAssistant< DIM, DataType >::phi(int s, int q, int var) const 
{
//  const size_t fe_ind = this->var_2_fe_type(var);
//  const size_t comp = this->var_2_comp(var);
//  return phi_[q][fe_ind][ this->fe_types_[ fe_ind ]->iv2ind(s, comp) ];
  return phi_[q][phi_offsets_for_var_[var] + s];
}

template < int DIM, class DataType >
inline DataType AssemblyAssistant< DIM, DataType >::Phi(int i, int q, int var) const 
{
  assert (var >= 0);
  assert (var < this->num_var_);
  assert (i >= 0);
  assert (i < this->num_dofs_total() );
  const size_t first_index = this->first_dof_for_var(var);
  
  if (i >= first_index && i < this->last_dof_for_var(var))
  {
    return this->phi(i - first_index, q, var);
  }
  return 0.;
}

template < int DIM, class DataType >
inline const Vec< DIM, DataType > & AssemblyAssistant< DIM, DataType >::grad_phi(int s, int q, int var) const 
{
//  const size_t fe_ind = this->var_2_fe_type(var);
//  const size_t comp = this->var_2_comp(var);
//  return grad_phi_[q][fe_ind][this->fe_types_[ fe_ind ]->iv2ind(s, comp)];
  return grad_phi_[q][phi_offsets_for_var_[var] + s];
}

template < int DIM, class DataType >
inline const Vec< DIM, DataType > & AssemblyAssistant< DIM, DataType >::grad_Phi(int i, int q, int var) const 
{
  assert (var >= 0);
  assert (var < this->num_var_);
  assert (i >= 0);
  assert (i < this->num_dofs_total() );
  const size_t first_index = this->first_dof_for_var(var);
  
  if (i >= first_index && i < this->last_dof_for_var(var))
  {
    return this->grad_phi(i - first_index, q, var);
  }
  return this->zero_vec_;
}

template < int DIM, class DataType >
inline const Mat< DIM, DIM, DataType > & AssemblyAssistant< DIM, DataType >::H_phi(int s, int q, int var) const 
{
  assert (H_phi_.size() > 0);
  return H_phi_[q][phi_offsets_for_var_[var] + s];
}

template < int DIM, class DataType >
inline DataType AssemblyAssistant< DIM, DataType >::L_phi(int s, int q, int var) const 
{
  assert (L_phi_.size() > 0);
  return L_phi_[q][phi_offsets_for_var_[var] + s];
}

template < int DIM, class DataType >
inline const Mat< DIM, DIM, DataType > & AssemblyAssistant< DIM, DataType >::H_Phi(int i, int q, int var) const 
{
  assert (var >= 0);
  assert (var < this->num_var_);
  assert (i >= 0);
  assert (i < this->num_dofs_total() );
  const size_t first_index = this->first_dof_for_var(var);
  
  if (i >= first_index && i < this->last_dof_for_var(var))
  {
    return this->H_phi(i - first_index, q, var);
  }
  return this->zero_mat_;
}

template < int DIM, class DataType >
inline DataType AssemblyAssistant< DIM, DataType >::L_Phi(int i, int q, int var) const 
{
  assert (var >= 0);
  assert (var < this->num_var_);
  assert (i >= 0);
  assert (i < this->num_dofs_total() );
  const size_t first_index = this->first_dof_for_var(var);
  
  if (i >= first_index && i < this->last_dof_for_var(var))
  {
    return this->L_phi(i - first_index, q, var);
  }
  return 0.;
}

template < int DIM, class DataType >
inline const Vec< DIM, DataType > & AssemblyAssistant< DIM, DataType >::x(int q) const 
{
  return x_[q];
}

template < int DIM, class DataType >
inline const Mat< DIM, DIM, DataType > & AssemblyAssistant< DIM, DataType >::J(int q) const 
{
  return J_[q];
}

template < int DIM, class DataType >
inline const Mat< DIM, DIM, DataType > & AssemblyAssistant< DIM, DataType >::JinvT(int q) const 
{
  return JinvT_[q];
}

template < int DIM, class DataType >
inline const std::vector< Mat< DIM, DIM, DataType > > & AssemblyAssistant< DIM, DataType >::H_F(int q) const 
{
  assert (H_.size() > 0);
  return H_[q];
}

template < int DIM, class DataType >
inline DataType AssemblyAssistant< DIM, DataType >::detJ(int q) const 
{
  return detJ_[q];
}

template < int DIM, class DataType >
inline const Mat< DIM, DIM - 1, DataType > & AssemblyAssistant< DIM, DataType >::Jf(int q) const 
{
  assert(facet_number_ > -1);
  return Jf_[q];
}

template < int DIM, class DataType >
inline DataType AssemblyAssistant< DIM, DataType >::ds(int q) const 
{
  assert(facet_number_ > -1);
  return ds_[q];
}

template < int DIM, class DataType >
inline const Vec< DIM, DataType > & AssemblyAssistant< DIM, DataType >::n(int q) const 
{
  assert(facet_number_ > -1);
  return mapped_n_[q];
}

template < int DIM, class DataType >
inline DataType AssemblyAssistant< DIM, DataType >::nitsche_regularization(int q) const 
{
  assert(facet_number_ > -1);
  // TODO: use q and specify this parameter somehow
  return nitsche_regularization_;
}

template < int DIM, class DataType >
inline const Vec< DIM, DataType > & AssemblyAssistant< DIM, DataType >::q_point(int q) const 
{
  return q_points_[q];
}

template < int DIM, class DataType >
inline int AssemblyAssistant< DIM, DataType >::num_quadrature_points() const 
{
  return q_points_.size();
}

template < int DIM, class DataType >
inline DataType AssemblyAssistant< DIM, DataType >::w(int q) const 
{
  return weights_[q];
}

template < int DIM, class DataType >
inline DataType AssemblyAssistant< DIM, DataType >::h() const 
{
  return h_;
}

template < int DIM, class DataType >
bool AssemblyAssistant< DIM, DataType >::need_cell_trafo_hessian (bool compute_hessians) const 
{
  bool ret = false;
  for (size_t i=0; i<fe_types_.size(); ++i)
  {
    if (fe_types_[i]->fe_trafo()->need_cell_trafo_hessian_for_gradients ())
    {
      return true;
    }
    if (compute_hessians)
    {
      if (fe_types_[i]->fe_trafo()->need_cell_trafo_hessian_for_hessians ())
      {
        return true;
      }
    }
  }
  return ret;
}

template < int DIM, class DataType >
bool AssemblyAssistant< DIM, DataType >::need_cell_trafo_grad_inv_detJ(bool compute_hessians) const 
{
  bool ret = false;
  for (size_t i=0; i<fe_types_.size(); ++i)
  {
    if (fe_types_[i]->fe_trafo()->need_cell_trafo_grad_inv_detJ_for_gradients())
    {
      return true;
    }
    if (compute_hessians)
    {
      if (fe_types_[i]->fe_trafo()->need_cell_trafo_grad_inv_detJ_for_hessians ())
      {
        return true;
      }
    }
  }
  return ret;
}

template < int DIM, class DataType >
template < class VectorType >
void AssemblyAssistant< DIM, DataType >::evaluate_fe_function(const VectorType &coefficients, 
                                                              int var,
                                                              FunctionValues< DataType > &function_values) const 
{
  // extract dof-values corresponding to local variable
  this->extract_dof_values(coefficients, var, this->local_coefficients_);

  // compute function values
  function_values.compute(
      phi_, EvalFiniteElementFunction< DIM, DataType >(var,
                                                       phi_offsets_for_var_, 
                                                       this->local_coefficients_) );
}

template < int DIM, class DataType >
template < class VectorType >
void AssemblyAssistant< DIM, DataType >::evaluate_fe_function_gradients( const VectorType &coefficients, 
                                                                         int var,
                                                                         FunctionValues< Vec< DIM, DataType > > &function_gradients) const 
{
  // extract dof-values corresponding to local variable
  this->extract_dof_values(coefficients, var, this->local_coefficients_);

  // compute function values
  function_gradients.compute(
    grad_phi_, EvalFiniteElementFunctionT< DIM, DataType, Vec<DIM, DataType> >(var,
                                                                               phi_offsets_for_var_, 
                                                                               this->local_coefficients_) );
}

template < int DIM, class DataType >
template < class VectorType >
void AssemblyAssistant< DIM, DataType >::evaluate_fe_function_hessians( const VectorType &coefficients, 
                                                                        int var,
                                                                        FunctionValues< Mat< DIM, DIM, DataType > > &function_hessians) const 
{
  assert (fe_H_phi_.size() > 0);

  // extract dof-values corresponding to local variable
  this->extract_dof_values(coefficients, var, this->local_coefficients_);

  // compute function values
  function_hessians.compute(
    H_phi_, EvalFiniteElementFunctionT< DIM, DataType, Mat<DIM, DIM, DataType> >(var,
                                                                                 phi_offsets_for_var_, 
                                                                                 this->local_coefficients_) );
}

template < int DIM, class DataType >
template < class VectorType >
void AssemblyAssistant< DIM, DataType >::evaluate_fe_function_hessians_and_laplace( const VectorType &coefficients, 
                                                                                    int var,
                                                                                    FunctionValues< Mat< DIM, DIM, DataType > > &function_hessians,
                                                                                    FunctionValues< DataType > &function_laplace) const 
{
  this->evaluate_fe_function_hessians(coefficients, var, function_hessians); 
  function_laplace.compute(function_hessians, EvalTrace<DIM, DataType>()); 
}

template < int DIM, class DataType >
void AssemblyAssistant< DIM, DataType >::initialize_for_element( const Element< DataType, DIM > &element,
                                                                 const Quadrature< DataType > &quadrature,
                                                                 bool compute_hessians) 
{  
  this->initialize(element, quadrature, -1, false, compute_hessians);
}

template < int DIM, class DataType >
void AssemblyAssistant< DIM, DataType >::initialize_for_facet( const Element< DataType, DIM > &element,
                                                               const Quadrature< DataType > &quadrature, 
                                                               int facet_number,
                                                               bool compute_hessians) 
{
  this->initialize(element, quadrature, facet_number, true, compute_hessians);
}

template < int DIM, class DataType >
void AssemblyAssistant< DIM, DataType >::initialize( const Element< DataType, DIM > &element,
                                                     const Quadrature< DataType > &quadrature,
                                                     int facet_number,
                                                     bool is_facet,
                                                     bool compute_hessians) 
{  
  const bool update_all = false;
  
  const bool first_call = (this->cell_index_ == -1);
  const bool changed_cell = (this->cell_index_ != element.cell_index());
  const bool changed_facet = (this->facet_number_ != facet_number);
  this->cell_index_ = element.cell_index();
  this->facet_number_ = facet_number;

  assert (this->cell_index_ >= 0);
  assert (!is_facet || this->facet_number_ >= 0);
  
  // update data structures for fe basis on reference cell
  const bool changed_fe_types = this->update_fe_types(element, update_all);
  assert (!first_call || changed_fe_types);
    
  // update quadrature on reference cell
  // TODO (Philipp G): for some reason, we have to force a quad update in case of DG assembly.
  // -> check whether this is really necessary or can be avoided
  const bool changed_quad_points = this->update_ref_quadrature(quadrature, is_facet || update_all);
  assert (!first_call || changed_quad_points);
  
  // evaluate fe basis functions on reference cell
  bool changed_ref_fe_val = false;
  if (changed_quad_points || changed_fe_types || update_all)
  {
    this->compute_ref_cell_values(compute_hessians);
    changed_ref_fe_val = true;
  }
  assert (!first_call || changed_ref_fe_val);

  // update cell transformation
  bool changed_cell_trafo = false;
  if (changed_cell || changed_quad_points || update_all)
  {
    changed_cell_trafo = this->update_cell_trafo(element, changed_quad_points || update_all, compute_hessians);
    // if false, new cell trafo differs from the old one by translation only 
    // and the quadrature points on the reference cells have not changed w.r.t. the previous call
  }
  assert (!first_call || changed_cell_trafo);
  
  // compute transformed values for new cell
  if (changed_ref_fe_val || changed_cell_trafo || update_all)
  {
    this->compute_trans_cell_values(compute_hessians);
  }
  
  // compute transformed values for current facet
  if (is_facet && (changed_facet || changed_fe_types || changed_cell_trafo || update_all))
  {
    // TODO: Compute this parameter in some way
    nitsche_regularization_ = 1.0;
    compute_facet_projection_matrix();
    compute_trans_facet_values();
    compute_facet_normal();
  }
  
  // update global dof indices
  this->global_dof_indices_.clear();
  element.get_dof_indices(this->global_dof_indices_);
  
  this->dof_factors_.clear();
  this->space_->dof().get_dof_factors_on_cell(element.cell_index(), this->dof_factors_);
}

template < int DIM, class DataType >
bool AssemblyAssistant< DIM, DataType >::update_fe_types( const Element< DataType, DIM > &element, 
                                                          bool force_update) 
{
  assert (element.get_space().fe_manager().is_initialized());
  
  const size_t num_fe = element.nb_fe();
  const size_t num_var = element.nb_var();
  
  const VectorSpace< DataType, DIM > *space = &element.get_space();
  const bool space_changed = !(space == this->space_);

  // fe_types_.empty() is true the first time function is called
  bool need_fe_update = force_update || fe_types_.empty() || space_changed;

  // check if element type changed
  if (!need_fe_update) 
  {
    for (size_t fe = 0; fe < num_fe; ++fe) 
    {
      const doffem::RefElement< DataType, DIM > *element_fe = element.get_fe(fe);

      assert(element_fe != 0);

      // check if element_fe_type already exists
      const int pos = find_fe_type(element_fe);
      if (pos < 0) 
      {
        need_fe_update = true;
        break;
      }
    }
  }

  if (!need_fe_update)
  {
    return false;
  }


  // clear and allocate data structures
  this->space_ = space;
  this->ref_cell_ = element.ref_cell();

  fe_types_.clear();
  fe_types_.reserve(num_fe);

  var_2_comp_.clear(); 
  var_2_fe_.clear();
  var_2_fe_type_.clear();
  fe_2_fe_type_.clear();    
  num_dofs_for_fe_.clear();   
  dof_offsets_for_fe_.clear();
  phi_offsets_for_var_.clear();
  //fe_offsets_for_comp_.clear();

  first_dof_for_var_.clear();
  last_dof_for_var_.clear();
    
  num_var_ = 0;
  num_dofs_ = 0;
  size_t var_offset = 0;
  size_t offset = 0;

  std::vector< size_t > fe_type2dof_offset;
  
  // loop over all elements present in FE space restricted to current cell
  for (size_t fe = 0; fe != num_fe; ++fe) 
  {
    doffem::RefElement< DataType, DIM > const * element_fe = element.get_fe(fe);
    assert(element_fe != 0);

    const size_t nb_comp = element_fe->nb_comp();
    const size_t dim = element_fe->dim();
      
    // avoid multiple evaluations of the same FE ansatz functions 
    // fe_types_ is list of unique RefElements  
    int fe_type = this->find_fe_type(element_fe);

    if (fe_type < 0) 
    { // this fe_type has not been added before
      fe_types_.push_back(element_fe);
      fe_type = fe_types_.size() - 1;
    } 

    num_dofs_for_fe_.push_back(element.nb_dof(fe));
    dof_offsets_for_fe_.push_back(offset);
    fe_2_fe_type_.push_back(fe_type);
      
    for (size_t comp = 0; comp < nb_comp; ++comp)
    {
      var_2_comp_.push_back(comp);
      var_2_fe_.push_back(fe);
      var_2_fe_type_.push_back(fe_type);
      first_dof_for_var_.push_back(num_dofs_);
      last_dof_for_var_.push_back(num_dofs_+dim);
        
      phi_offsets_for_var_.push_back(var_offset);
      var_offset += dim;
    }
    num_var_ += nb_comp;

    offset += element.nb_dof(fe);
    num_dofs_ += dim;
  }
  phi_offsets_for_var_.push_back(var_offset);

  this->initialized_fe_types_ = true;
  return true;
}

template < int DIM, class DataType >
bool AssemblyAssistant< DIM, DataType >::update_ref_quadrature( const Quadrature< DataType > &new_quadrature, 
                                                                bool force_update) 
{ 
  // check if quadrature changed
  bool need_quadrature_update =  force_update 
                              || (quadrature_.size() == 0) 
                              || (new_quadrature.size() != quadrature_.size())
                              || (new_quadrature.name() != quadrature_.name());

  if (!need_quadrature_update) 
  {
    return false;
  }
  
  // copy quadrature
  quadrature_ = new_quadrature;
    
  this->q_points_prev_.clear();
  this->q_points_prev_.insert(this->q_points_prev_.begin(), this->q_points_.begin(),this->q_points_.end());
  
  const size_t num_q = quadrature_.size();
  q_points_.resize(num_q);
  weights_.resize(num_q);
  for (size_t q = 0; q != num_q; ++q) 
  {
    for (size_t c = 0; c != DIM; ++c) 
    {
      switch (c) 
      {
        case 0:
          q_points_[q][0] = quadrature_.x(q);
          break;
        case 1:
          q_points_[q][1] = quadrature_.y(q);
          break;
        case 2:
          q_points_[q][2] = quadrature_.z(q);
          break;
        default:
          assert(false);
          break;
      }
    }
    weights_[q] = quadrature_.w(q);
  }

  if (this->q_points_prev_.size() == 0)
  {
    // first call to this routine
    return true;
  }
  return !(vectors_are_equal<DataType, DIM>(this->q_points_, this->q_points_prev_, 1e-10));
}

template < int DIM, class DataType >
bool AssemblyAssistant< DIM, DataType >::update_cell_trafo(const Element< DataType, DIM > &element,
                                                           bool changed_quad,
                                                           bool compute_hessians) 
{
  bool differs_by_translation_only = false;
  
  // check whether cell transformation differs from previous one only by translation
  if (this->cell_transform_)
  {
    // this is not the first time, this routine is called
    assert (cell_transform_.get() != nullptr);
    differs_by_translation_only = (this->cell_transform_->differs_by_translation_from(element.get_cell_transformation()));
  }
  this->cell_transform_ = element.get_cell_transformation();
  
  assert (q_points_.size() == quadrature_.size());
  
  // compute quadrature points on physical cell 
  x_.compute(q_points_, EvalPhysicalPoint< DIM, DataType >(*cell_transform_));
  
  if (differs_by_translation_only && (!changed_quad))
  {
    return false;
  }
  
  // compute jacobians on physical cell 
  J_.compute(q_points_, EvalPhysicalJacobian< DIM, DataType >(*cell_transform_));

  // compute determinant of jacobians 
  detJ_.compute(J_, EvalDeterminant< DIM, DataType >());

  // compute gradient of inverse determinant of jacobians 
  if (this->need_cell_trafo_grad_inv_detJ(compute_hessians))
  {
    grad_inv_detJ_.compute(q_points_, EvalGradInvDeterminantCellTrafo< DIM, DataType >(*cell_transform_));
  }
  else
  {
    grad_inv_detJ_.zeros(q_points_.size());
  }
  
  // compute inverse of jacobians
  Jinv_.compute(J_, EvalInverse< DIM, DataType >());
  
  // compute transpose of inverse jacobians
  JinvT_.compute(Jinv_, EvalTranspose< DIM, DataType >());

  if (this->need_cell_trafo_hessian(compute_hessians))
  {
    // compute hessians on physical cell 
    H_.compute(q_points_, EvalCellTransformationHessian< DIM, DataType >(*cell_transform_));
  }
  else
  {
    H_.zeros(q_points_.size());
  }
  
  // compute mesh parameter h
  h_ = std::pow(static_cast< DataType >(std::abs(detJ(0))),
                  static_cast< DataType >(1. / static_cast< DataType >(DIM)));

  return true;
}

template < int DIM, class DataType >
void AssemblyAssistant< DIM, DataType >::compute_ref_cell_values(bool compute_hessians) 
{
  assert (q_points_.size() > 0);
  assert(q_points_.size() == quadrature_.size());
  assert(!fe_types_.empty());
  assert (initialized_fe_types_);
  
  // compute shape function values 
  phi_hat_.compute(q_points_, EvalShapeFunctions< DIM, DataType >(fe_types_));

  // compute shape function gradients 
  grad_phi_hat_.compute(q_points_, EvalShapeFunctionGradients< DIM, DataType >(fe_types_));

  if (compute_hessians)
  {
    // compute shape function hessians 
    H_phi_hat_.compute(q_points_, EvalShapeFunctionHessians< DIM, DataType >(fe_types_));
  }
}

template < int DIM, class DataType >
void AssemblyAssistant< DIM, DataType >::compute_trans_cell_values(bool compute_hessians) 
{
  assert (q_points_.size() == quadrature_.size());
  assert (cell_transform_ != 0);
  assert (detJ_.size() == J_.size());
  assert (detJ_.size() == JinvT_.size());
  assert (detJ_.size() == Jinv_.size());
  assert (detJ_.size() == grad_inv_detJ_.size());
  assert (detJ_.size() == H_.size());
  assert (detJ_.size() == phi_hat_.size());
  assert (detJ_.size() == grad_phi_hat_.size());
      
  // compute shape function values 
  fe_phi_.compute(detJ_, J_, JinvT_,
                    EvalMappedShapeFunctions< DIM, DataType >(fe_types_, phi_hat_));

  phi_.compute (detJ_, ReorderMappedShapeT< DIM, DataType, DataType >
                        (fe_types_, var_2_fe_type_, var_2_comp_, fe_phi_, phi_offsets_for_var_));

  // compute JinvT * grad_phi_hat 
  fe_grad_phi_.compute (detJ_, grad_inv_detJ_, J_, Jinv_, JinvT_, H_, 
                          EvalMappedShapeFunctionGradients< DIM, DataType >(fe_types_, fe_phi_, phi_hat_, grad_phi_hat_));

  grad_phi_.compute (detJ_, ReorderMappedShapeT< DIM, DataType, Vec<DIM, DataType> >
                                (fe_types_, var_2_fe_type_, var_2_comp_, fe_grad_phi_, phi_offsets_for_var_));
  
  if (compute_hessians)
  {
    fe_H_phi_.compute (JinvT_, H_, 
                         EvalMappedShapeFunctionHessians< DIM, DataType >( fe_types_, fe_grad_phi_, grad_phi_hat_, H_phi_hat_));
  
    H_phi_.compute (detJ_, ReorderMappedShapeT< DIM, DataType, Mat<DIM, DIM, DataType> >
                              (fe_types_, var_2_fe_type_, var_2_comp_, fe_H_phi_, phi_offsets_for_var_));

    L_phi_.compute_wo_arg(H_phi_.size(), EvalLaplacian<DIM, DataType>(num_var_, H_phi_));
  }
}

template < int DIM, class DataType >
void AssemblyAssistant< DIM, DataType >::compute_facet_projection_matrix() 
{
  assert (this->ref_cell_);
  this->ref_cell_->compute_facet_projection_matrix(this->facet_number_, this->facet_proj_);
}

template < int DIM, class DataType >
void AssemblyAssistant< DIM, DataType >::compute_trans_facet_values() 
{
  assert(q_points_.size() == quadrature_.size());
  assert(cell_transform_ );

  Jf_.compute(J_, EvalRightMatrixMult< DIM, DIM, DIM - 1, DataType >(facet_proj_));
  ds_.compute(Jf_, EvalSurfaceElement< DIM, DataType >());
}

template < int DIM, class DataType >
void AssemblyAssistant< DIM, DataType >::compute_facet_normal() 
{   
  assert (this->ref_cell_);
  this->ref_cell_->compute_facet_normal(this->facet_number_, this->n_);
 
  // precompute the normals on the quadrature points.
  mapped_n_.compute(JinvT_, EvalMappedNormal< DIM, DataType >(n_));
}

template < int DIM, class DataType >
int AssemblyAssistant< DIM, DataType >::find_fe_type(const doffem::RefElement< DataType, DIM > *fe_type) const 
{
  assert(fe_type != 0);
  int pos = 0;
  bool fe_type_found = false;

  //typedef typename std::vector< doffem::RefElement< DataType, DIM > const * >::const_iterator Iterator;
  for (size_t i=0, e = this->fe_types_.size(); i!=e; ++i) 
  {
    if ((*this->fe_types_[i]) == (*fe_type)) 
    {
      fe_type_found = true;
      break;
    }
    ++pos;
  }

  if (fe_type_found) 
  {
    return pos;
  }
  return -1;
}

// TODO: optimize use of std::vector
template < int DIM, class DataType >
template < class VectorType >
void AssemblyAssistant< DIM, DataType >::extract_dof_values( const VectorType &coefficients, 
                                                             size_t var,
                                                             std::vector<DataType>& local_coeff) const 
{
  assert (var < this->var_2_fe_.size());
  const size_t fe_ind = this->var_2_fe (var);

  // extract dof-values corresponding to local finite element

  // find global dof indices corresponding to this element
  const size_t num_fe_dofs = num_dofs_for_fe(fe_ind);
  
  local_coeff.clear();
  local_coeff.resize(num_fe_dofs, 0.);
  
  this->tmp_indices_.clear();
  this->tmp_indices_.resize(num_fe_dofs, -1);
  
  this->tmp_factors_.clear();
  this->tmp_factors_.resize(num_fe_dofs, 1.);
    
  for (size_t i = 0; i != num_fe_dofs; ++i) 
  {
    this->tmp_indices_[i] = this->global_dof_indices_[this->dof_index_for_fe(i, fe_ind)];  
    this->tmp_factors_[i] = this->dof_factors_[this->dof_index_for_fe(i, fe_ind)];
  }
  coefficients.GetValues(vec2ptr(this->tmp_indices_), num_fe_dofs, vec2ptr(local_coeff));
  
  for (size_t i = 0; i != num_fe_dofs; ++i) 
  {
    local_coeff[i] *= this->tmp_factors_[i]; 
  }
}

} // namespace hiflow

#endif /* _ASSEMBLY_ASSISTANT_H_ */
