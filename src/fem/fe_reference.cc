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

#include "fem/fe_reference.h"
#include "dof/dof_impl/dof_container.h"
#include "fem/ansatz/ansatz_space.h"
#include "linear_algebra/seq_dense_matrix.h"
#include <iomanip>

namespace hiflow {
namespace doffem {

template < class DataType, int DIM >
RefElement< DataType, DIM >::RefElement()
    : 
    space_(nullptr), 
    dofs_(nullptr),
    instance_id_(-1), 
    type_(FE_TYPE_NOT_SET), 
    fe_trafo_(nullptr),
    is_nodal_basis_(false),
    V_inv_(nullptr),
    V_(nullptr) 
{}

template < class DataType, int DIM > 
RefElement< DataType, DIM >::~RefElement() 
{
  if (this->V_inv_ != nullptr)
  {
    delete this->V_inv_;
  }
  if (this->V_ != nullptr)
  {
    delete this->V_;
  }
}

template < class DataType, int DIM >
size_t RefElement< DataType, DIM >::nb_subentity(int tdim) const 
{
  assert(this->dofs_ != nullptr);
  return this->dofs_->nb_subentity(tdim);
}

template < class DataType, int DIM > 
size_t RefElement< DataType, DIM >::nb_dof_on_cell() const 
{
  assert(this->dofs_ != nullptr);
  return this->dofs_->nb_dof_on_cell();
}

template < class DataType, int DIM >
size_t RefElement< DataType, DIM >::nb_dof_on_subentity(int tdim, int index) const 
{
  assert(this->dofs_ != nullptr);
  return this->dofs_->nb_dof_on_subentity(tdim, index);
}

template < class DataType, int DIM >
std::vector< DofID > const & RefElement< DataType, DIM >::get_dof_on_subentity(int tdim, int index) const 
{
  assert(this->dofs_ != nullptr);
  return this->dofs_->get_dof_on_subentity(tdim, index);
}

template < class DataType, int DIM >
void RefElement< DataType, DIM >::init (AnsatzSpace<DataType, DIM> const * space, 
                                        DofContainer<DataType, DIM> const * dofs,
                                        FETransformation<DataType, DIM> const * fe_trafo,
                                        bool is_nodal_basis,
                                        FEType fe_type)
{
  assert (space != nullptr);
  assert (dofs != nullptr);
  assert (fe_trafo != nullptr);
  assert (fe_type != FE_TYPE_NOT_SET);
  
  // set some parameters
  this->ref_cell_ = space->ref_cell();
  this->dim_ = space->dim();
  this->tdim_ = space->tdim();
  this->max_deg_ = space->max_deg();
  this->nb_comp_ = space->nb_comp();
  this->weight_size_ = space->weight_size();
  this->ref_cell_ = space->ref_cell();
  this->type_ = fe_type;
    
/*
  this->comp_weight_size_.resize(this->weight_size);
  for (size_t i=0; i<this->weight_size_; ++i)
  {
    this->comp_weight_size_[i] = space->get_weight_size(i);
  }
*/

  this->fe_trafo_ = fe_trafo;
  this->is_nodal_basis_ = is_nodal_basis;
  this->dofs_ = dofs;
  this->space_ = space;

  // check that dofcontainer and ansatz space are compatible
  assert (dofs->nb_dof_on_cell() == this->dim_);
  assert (dofs->ref_cell_type() == space->ref_cell_type());
    
  // compute transformation matrix V_inv:
  // phi_i (x) = sum_k V_inv(i,k) * psi_k(x)
  // where psi_k is the k-th basis function defined in space_
  // and phi_i is the i-th basis function such that
  // dof_j (phi_i) = delta_{i,j} holds
    
  if (!is_nodal_basis_)
  {
    this->compute_basis_transformation_matrix();
  }

  this->name_ = space->name() + "_" + dofs->name();
    
  this->init_status_ = true;
}

template < class DataType, int DIM >
bool RefElement< DataType, DIM >::operator==(const RefElement< DataType, DIM > &fe_slave) const 
{
  if (this->type_ == FE_TYPE_NOT_SET || fe_slave.type_ == FE_TYPE_NOT_SET) 
  {
    std::cout << this->type_ << " =? " << fe_slave.type_ << " " << FE_TYPE_NOT_SET << std::endl;
    assert(0);
  }
  assert (this->space_ != nullptr);
  assert (fe_slave.space_ != nullptr);
  assert (this->dofs_ != nullptr);
  assert (fe_slave.dofs_ != nullptr);
  
  return this->type_ == fe_slave.type_ 
      && (*(this->space_) == *(fe_slave.space_))
      && (*(this->dofs_) == *(fe_slave.dofs_));
}

template < class DataType, int DIM >
bool RefElement< DataType, DIM >::operator<(const RefElement< DataType, DIM > &fe_slave) const 
{
  if (this->type_ == FE_TYPE_NOT_SET || fe_slave.type_ == FE_TYPE_NOT_SET) {
    assert(0);
  }
  assert (this->space_ != nullptr);
  assert (fe_slave.space_ != nullptr);
  assert (this->dofs_ != nullptr);
  assert (fe_slave.dofs_ != nullptr);


  if (this->type_ < fe_slave.type_) 
  {
    return true;
  } 
  else if (this->type_ == fe_slave.type_) 
  {
    if (*(this->space_) < *(fe_slave.space_))
    {
      return true;
    }
    else if (*(this->space_) == *(fe_slave.space_))
    {
      if (*(this->dofs_) < *(fe_slave.dofs_))
      {
         return true;
      }
    }
  }
  return false;
}

template < class DataType, int DIM >
void RefElement< DataType, DIM >::compute_degree_hash ( ) const
{
  assert (this->space_ != nullptr);
//this->space_->compute_degree_hash();
  this->deg_hash_ = this->space_->deg_hash();
}

template < class DataType, int DIM >
void RefElement< DataType, DIM >::compute_basis_transformation_matrix ( ) 
{
  assert (this->space_ != nullptr);
  assert (this->dofs_->nb_dof_on_cell() == this->space_->dim());
  assert (this->dim() == this->space_->dim());
  
  if (this->V_inv_ != nullptr)
  {
    delete this->V_inv_;
  }
  if (this->V_ != nullptr)
  {
    delete this->V_;
  }

  this->V_ = new la::SeqDenseMatrix<DataType>();
  this->V_inv_ = new la::SeqDenseMatrix<DataType>();

  const size_t n = this->dim();
  this->V_->Resize(n,n);
  
  // compute basis trafo matrix: V_{j,i} = dof_i( space->basis_j )
  std::vector< DofID > all_dofs(n, 0);
  for (size_t i = 0; i<n; ++i)
  {
    all_dofs[i] = i;
  }
  
  assert (this->dofs_ != nullptr);

  std::vector< std::vector<DataType> >dof_values;
  this->dofs_->evaluate(this->space_, all_dofs, dof_values);
    
  assert (dof_values.size() == all_dofs.size());
  
  for (size_t i = 0; i<n; ++i)
  {
    assert (dof_values[i].size() == n);
    for (size_t j=0; j<n; ++j)
    {
      this->V_->operator()(j,i) = dof_values[i][j];
    }
  }

/*
  for (size_t i = 0; i<n; ++i)
  {
    for (size_t j=0; j<n; ++j)
    {
      DataType val = this->V_->operator()(i,j);
      if (std::abs(val) < 1e-14)
      {
        val = 0.;
      }
      std::cout << std::fixed << std::setw( 7 ) << std::setprecision( 4 ) << val << " ";
    }
    std::cout << std::endl;
  }
*/
  // compute inverse of V
  bool success = this->V_->Factorize();
  assert (success);
  
  this->V_inv_->Resize(n,n);
  for (size_t i=0; i<n; ++i)
  {
    std::vector< DataType > b (n, 0.);
    std::vector< DataType > x (n, 0.);
    b[i] = 1.;
  
    // V * x = b = e_i
    this->V_->ForwardBackward(b, x);

    // V_inv [:,i] = x
    for (size_t j=0; j<n; ++j)
    {
      this->V_inv_->operator()(j,i) = x[j];
    }
  }
}

template < class DataType, int DIM >
void RefElement< DataType, DIM >::apply_transformation (const std::vector< DataType > &weight_in, 
                                                        std::vector< DataType > &weight_out) const
{
  assert (weight_in.size() == weight_out.size());
  assert (this->init_status_);
  
  const size_t n = this->dim();
  
  // phi_i (x) = sum_k V_inv(i,k) * psi_k(x)
  
  // loop over components
  for (size_t v=0, e = this->nb_comp(); v!=e; ++v)
  {
    // loop over new basis functions
    for (size_t i = 0; i < n; ++i)
    { 
      weight_out[this->iv2ind(i,v)] = 0.;
      
      // loop over old basis functions
      for (size_t k = 0; k < n; ++k)
      {
        weight_out[this->iv2ind(i,v)] += this->V_inv_->operator()(i,k) * weight_in[this->iv2ind(k,v)];
      }
    }
  }
}

template < class DataType, int DIM >
void RefElement< DataType, DIM >::evaluate_and_transform(BasisEvalFunction fun,  
                                                         const Coord &pt, 
                                                         std::vector< DataType > &weight) const
{
  assert (this->dim() * this->nb_comp() == weight.size());
  if (!this->is_nodal_basis_)
  {
    // AnsatzSpace basis does not satisfy nodal property -> transformation needed
    this->weight_psi_.clear();
    this->weight_psi_.resize(this->dim() * this->nb_comp(), 0.);
    fun(pt, this->weight_psi_);
    this->apply_transformation(this->weight_psi_, weight);
  }
  else
  {
    // AnsatzSpace basis is already nodal by definition -> no transformation needed
    fun(pt, weight);
  }
}

template < class DataType, int DIM >
void RefElement< DataType, DIM >::N(const Coord &pt, std::vector< DataType > &weight) const 
{
  BasisEvalFunction fun = boost::bind ( &AnsatzSpace< DataType, DIM >::N, this->space_, _1, _2 );
  this->evaluate_and_transform(fun, pt, weight);
//std::cout << "[RefElement] " << pt[0] << ", " << pt[1] << ": " << weight[0] << " " << weight[1] << "..." << std::endl;
}

template < class DataType, int DIM >
void RefElement< DataType, DIM >::N_x(const Coord &pt, std::vector< DataType > &weight) const 
{
  BasisEvalFunction fun = boost::bind ( &AnsatzSpace< DataType, DIM >::N_x, this->space_, _1, _2 );
  this->evaluate_and_transform(fun, pt, weight);
}

template < class DataType, int DIM >
void RefElement< DataType, DIM >::N_y(const Coord &pt, std::vector< DataType > &weight) const 
{
  BasisEvalFunction fun = boost::bind ( &AnsatzSpace< DataType, DIM >::N_y, this->space_, _1, _2 );
  this->evaluate_and_transform(fun, pt, weight);
}

template < class DataType, int DIM >
void RefElement< DataType, DIM >::N_z(const Coord &pt, std::vector< DataType > &weight) const 
{
  BasisEvalFunction fun = boost::bind ( &AnsatzSpace< DataType, DIM >::N_z, this->space_, _1, _2 );
  this->evaluate_and_transform(fun, pt, weight);
}

template < class DataType, int DIM >
void RefElement< DataType, DIM >::N_xx(const Coord &pt, std::vector< DataType > &weight) const 
{
  BasisEvalFunction fun = boost::bind ( &AnsatzSpace< DataType, DIM >::N_xx, this->space_, _1, _2 );
  this->evaluate_and_transform(fun, pt, weight);
}

template < class DataType, int DIM >
void RefElement< DataType, DIM >::N_xy(const Coord &pt, std::vector< DataType > &weight) const 
{
  BasisEvalFunction fun = boost::bind ( &AnsatzSpace< DataType, DIM >::N_xy, this->space_, _1, _2 );
  this->evaluate_and_transform(fun, pt, weight);
}

template < class DataType, int DIM >
void RefElement< DataType, DIM >::N_xz(const Coord &pt, std::vector< DataType > &weight) const 
{
  BasisEvalFunction fun = boost::bind ( &AnsatzSpace< DataType, DIM >::N_xz, this->space_, _1, _2 );
  this->evaluate_and_transform(fun, pt, weight);
}

template < class DataType, int DIM >
void RefElement< DataType, DIM >::N_yy(const Coord &pt, std::vector< DataType > &weight) const 
{
  BasisEvalFunction fun = boost::bind ( &AnsatzSpace< DataType, DIM >::N_yy, this->space_, _1, _2 );
  this->evaluate_and_transform(fun, pt, weight);
}

template < class DataType, int DIM >
void RefElement< DataType, DIM >::N_yz(const Coord &pt, std::vector< DataType > &weight) const 
{
  BasisEvalFunction fun = boost::bind ( &AnsatzSpace< DataType, DIM >::N_yz, this->space_, _1, _2 );
  this->evaluate_and_transform(fun, pt, weight);
}

template < class DataType, int DIM >
void RefElement< DataType, DIM >::N_zz(const Coord &pt, std::vector< DataType > &weight) const 
{
  BasisEvalFunction fun = boost::bind ( &AnsatzSpace< DataType, DIM >::N_zz, this->space_, _1, _2 );
  this->evaluate_and_transform(fun, pt, weight);
}


template class RefElement< float, 3 >;
template class RefElement< float, 2 >;
template class RefElement< float, 1 >;

template class RefElement< double, 3 >;
template class RefElement< double, 2 >;
template class RefElement< double, 1 >;

} // namespace doffem
} // namespace hiflow
