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

#include "space/element.h"

#include "dof/dof_partition.h"
#include "fem/cell_trafo/cell_transformation.h"
#include "fem/fe_reference.h"
#include "fem/fe_manager.h"
#include "fem/fe_transformation.h"
#include "space/vector_space.h"
#include "mesh/entity.h"
#include "mesh/iterator.h"

#include <boost/bind/bind.hpp> 

namespace hiflow {

using namespace boost::placeholders;

template <class DataType, int DIM>
Element<DataType, DIM>::Element(const VectorSpace< DataType, DIM > &space, int cell_index)
  : space_(&space),
    cell_(space.mesh().get_entity(space.tdim(), cell_index)) 
{
  this->dim_ = 0;
  this->nb_comp_ = 0;
  this->weight_size_ = 0; 
  this->mapping_eval_ = false;
  
  const size_t nb_fe = this->space_->fe_manager().nb_fe();
  this->weight_offsets_.resize(nb_fe, 0);
  this->dim_offsets_.resize(nb_fe, 0);

  size_t weight_offset = 0;
  size_t dim_offset = 0;
  for (size_t i=0; i<nb_fe; ++i)
  { 
    const doffem::RefElement< DataType, DIM > * cur_fe = this->space_->fe_manager().get_fe(this->cell_.index(), i);
      
    this->weight_offsets_[i] = weight_offset;
    this->dim_offsets_[i] = dim_offset;
    this->nb_comp_ += cur_fe->nb_comp();

    weight_offset += cur_fe->weight_size();
    dim_offset += cur_fe->dim();
  }
  this->weight_size_ = weight_offset;
  this->dim_ = dim_offset;
}

template <class DataType, int DIM>
std::shared_ptr< doffem::CellTransformation< DataType, DIM > > Element<DataType, DIM>::get_cell_transformation(/*size_t fe_ind*/) const 
{
  return this->space_->fe_manager().get_cell_transformation(this->cell_.index()/*, fe_ind*/);
}

template <class DataType, int DIM>
const doffem::RefElement< DataType, DIM > *Element<DataType, DIM>::get_fe(size_t fe_ind) const  
{
  return this->space_->fe_manager().get_fe(this->cell_.index(), fe_ind);
}

template <class DataType, int DIM>
const doffem::RefElement< DataType, DIM > *Element<DataType, DIM>::get_fe_for_var(size_t var) const 
{
  return this->space_->fe_manager().get_fe_for_var(this->cell_.index(), var);
}

template <class DataType, int DIM>
void Element<DataType, DIM>::get_dof_indices(size_t fe_ind, std::vector< int > &indices) const 
{
  this->space_->get_dof_indices(fe_ind, this->cell_.index(), &indices);
}

template <class DataType, int DIM>
void Element<DataType, DIM>::get_dof_indices(std::vector< int > &indices) const 
{
  this->space_->get_dof_indices(this->cell_.index(), &indices);
}

template <class DataType, int DIM>
void Element<DataType, DIM>::get_dof_indices_on_subentity(size_t fe_ind, int tdim, int sindex, std::vector< int > &indices) const 
{
  this->space_->get_dof_indices_on_subentity(fe_ind, this->cell_.index(), tdim, sindex, &indices);
}

template <class DataType, int DIM>
void Element<DataType, DIM>::get_dof_indices_on_subentity(int tdim, int sindex,
                                                          std::vector< int > &indices) const 
{
  this->space_->get_dof_indices_on_subentity(this->cell_.index(), tdim, sindex, &indices);
}

template <class DataType, int DIM>
size_t Element<DataType, DIM>::nb_var() const 
{
  return this->space_->fe_manager().nb_var();
}

template <class DataType, int DIM>
size_t Element<DataType, DIM>::nb_fe() const 
{
  return this->space_->fe_manager().nb_fe();
}

template <class DataType, int DIM>
size_t Element<DataType, DIM>::dof_offset(size_t fe_ind) const 
{
  return this->dim_offsets_[fe_ind];
}

template <class DataType, int DIM>
size_t Element<DataType, DIM>::var_2_fe (size_t var) const
{
  return this->space_->fe_manager().var_2_fe(var);
}

template <class DataType, int DIM>
size_t Element<DataType, DIM>::var_2_comp (size_t var) const
{
  return this->space_->fe_manager().var_2_comp(var);
}

template <class DataType, int DIM>
std::vector<size_t> Element<DataType, DIM>::fe_2_var (size_t fe_ind) const
{
  return this->space_->fe_manager().fe_2_var(fe_ind);
}

template <class DataType, int DIM>
doffem::ConstRefCellPtr<DataType, DIM> Element<DataType, DIM>::ref_cell() const 
{
  return this->space_->fe_manager().get_fe(this->cell_.index(), 0)->ref_cell();
}

template <class DataType, int DIM>
mesh::EntityNumber Element<DataType, DIM>::cell_index() const 
{ 
  return this->cell_.index(); 
}

template <class DataType, int DIM>
size_t Element<DataType, DIM>::nb_dof(size_t fe_ind) const 
{
  const doffem::RefElement< DataType, DIM > *fe_type = this->get_fe(fe_ind);
  assert(fe_type != 0);
  return fe_type->nb_dof_on_cell();
}

template <class DataType, int DIM>
size_t Element<DataType, DIM>::nb_comp (size_t fe_ind) const
{
  return this->space_->fe_manager().get_fe(this->cell_.index(), fe_ind)->nb_comp();
}
 
template<class DataType, int DIM>
bool Element<DataType, DIM>::is_boundary() const 
{
  const mesh::TDim facet_dim = this->space_->mesh().tdim() - 1;
  for (mesh::IncidentEntityIterator it = this->cell_.begin_incident(facet_dim);
       it != this->cell_.end_incident(facet_dim); ++it) 
  {
    bool is_boundary = this->space_->mesh().is_boundary_facet(it->index());
    if (is_boundary) 
    {
      return true;
    }
  }
  return false;
}


template<class DataType, int DIM>
bool Element<DataType, DIM>::is_local() const 
{
  const mesh::TDim tdim = this->space_->mesh().tdim();
  if (this->space_->mesh().has_attribute("_sub_domain_", tdim)) {
    int cell_sub_subdomain;
    this->cell_.get("_sub_domain_", &cell_sub_subdomain);
    return cell_sub_subdomain == this->space_->dof().my_subdom();
  }
  // assume it is true if we have no other information
  return true;
}

template<class DataType, int DIM>
std::vector< int > Element<DataType, DIM>::boundary_facet_numbers() const 
{
  const mesh::TDim facet_dim = this->space_->mesh().tdim() - 1;
  std::vector< int > facet_numbers(0);
  if (is_boundary()) {
    facet_numbers.reserve(6);
    int i = 0;
    for (mesh::IncidentEntityIterator it = this->cell_.begin_incident(facet_dim);
         it != this->cell_.end_incident(facet_dim); ++it) 
    {
      if (this->space_->mesh().is_boundary_facet(it->index())) 
      {
        facet_numbers.push_back(i);
      }
      ++i;
    }
  }
  return facet_numbers;
}

template<class DataType, int DIM>
size_t Element<DataType, DIM>::iv2ind (size_t i, size_t var) const
{
  const size_t fe_ind = this->var_2_fe(var);

#ifndef NDEBUG
  assert (this->var_2_comp(var) < this->space_->fe_manager().get_fe(this->cell_.index(), fe_ind)->nb_comp());
  assert (i < this->space_->fe_manager().get_fe(this->cell_.index(), fe_ind)->dim());
#endif

  return this->weight_offsets_[fe_ind] 
         + this->space_->fe_manager().get_fe(this->cell_.index(), fe_ind)->iv2ind(i, this->var_2_comp(var));
}

template<class DataType, int DIM>
void Element<DataType, DIM>::N (const Coord &ref_pt, std::vector< DataType > &weights) const
{ 
  const size_t nb_fe = this->nb_fe();
  size_t offset = 0;

  // loop over fe types
  for (size_t i=0; i<nb_fe; ++i)
  {
    doffem::RefElement<DataType, DIM> const * ref_fe = this->get_fe(i);
    const size_t cur_weight_size = ref_fe->weight_size();
    
    std::vector< DataType> cur_weights (cur_weight_size, 0.);
    this->N_fe (ref_pt, i, cur_weights);

    assert (weights.size() > offset + cur_weight_size);
    for (size_t l=0; l<cur_weight_size; ++l)
    {
      weights[offset+l] = cur_weights[l];
    }
    offset += cur_weight_size;
  }
}

template<class DataType, int DIM>
void Element<DataType, DIM>::grad_N (const Coord &ref_pt, std::vector< Vec<DIM,DataType> > &gradients) const
{ 
  const size_t nb_fe = this->nb_fe();
  size_t offset = 0;

  // loop over fe types
  for (size_t i=0; i<nb_fe; ++i)
  {
    doffem::RefElement<DataType, DIM> const * ref_fe = this->get_fe(i);
    const size_t cur_weight_size = ref_fe->weight_size();
    
    std::vector< Vec<DIM,DataType> > cur_grads (cur_weight_size);
    this->grad_N_fe (ref_pt, i, cur_grads);

    assert (gradients.size() > offset + cur_weight_size);
    for (size_t l=0; l<cur_weight_size; ++l)
    {
      gradients[offset+l] = cur_grads[l];
    }
    offset += cur_weight_size;
  }
}

template<class DataType, int DIM>
void Element<DataType, DIM>::hessian_N (const Coord &ref_pt, std::vector< Mat<DIM, DIM, DataType> > &hessians) const
{ 
  const size_t nb_fe = this->nb_fe();
  size_t offset = 0;

  // loop over fe types
  for (size_t i=0; i<nb_fe; ++i)
  {
    doffem::RefElement<DataType, DIM> const * ref_fe = this->get_fe(i);
    const size_t cur_weight_size = ref_fe->weight_size();
    
    std::vector< Mat<DIM, DIM, DataType> > cur_hessians (cur_weight_size);
    this->hessian_N_fe (ref_pt, i, cur_hessians);

    assert (hessians.size() > offset + cur_weight_size);
    for (size_t l=0; l<cur_weight_size; ++l)
    {
      hessians[offset+l] = cur_hessians[l];
    }
    offset += cur_weight_size;
  }
}

template<class DataType, int DIM>
void Element<DataType, DIM>::N_var (const Coord &ref_pt, size_t var, std::vector< DataType > &weights) const
{ 
  assert (var < this->nb_var());
  const size_t fe_ind = this->var_2_fe(var);
  const size_t comp = this->var_2_comp(var);

  doffem::RefElement<DataType, DIM> const * ref_fe = this->get_fe(fe_ind);
  const size_t fe_weight_size = ref_fe->weight_size();
  const size_t dim = ref_fe->dim();

  std::vector< DataType> fe_weights (fe_weight_size, 0.);
  this->N_fe (ref_pt, fe_ind, fe_weights);

  assert (weights.size() == dim);
  for (size_t i=0; i<dim; ++i)
  {
    weights[i] = fe_weights[ref_fe->iv2ind(i,comp)];
  }
}

template<class DataType, int DIM>
void Element<DataType, DIM>::grad_N_var (const Coord &ref_pt, size_t var, std::vector< Vec<DIM,DataType> > &gradients) const
{ 
  assert (var < this->nb_var());
  const size_t fe_ind = this->var_2_fe(var);
  const size_t comp = this->var_2_comp(var);

  doffem::RefElement<DataType, DIM> const * ref_fe = this->get_fe(fe_ind);
  const size_t fe_weight_size = ref_fe->weight_size();
  const size_t dim = ref_fe->dim();

  std::vector< Vec<DIM,DataType> > fe_weights (fe_weight_size);
  this->grad_N_fe (ref_pt, fe_ind, fe_weights);

  assert (gradients.size() == dim);
  for (size_t i=0; i<dim; ++i)
  {
    gradients[i] = fe_weights[ref_fe->iv2ind(i,comp)];
  }
}

template<class DataType, int DIM>
void Element<DataType, DIM>::hessian_N_var (const Coord &ref_pt, size_t var, std::vector< Mat<DIM,DIM,DataType> > &hessians) const
{ 
  assert (var < this->nb_var());
  const size_t fe_ind = this->var_2_fe(var);
  const size_t comp = this->var_2_comp(var);

  doffem::RefElement<DataType, DIM> const * ref_fe = this->get_fe(fe_ind);
  const size_t fe_weight_size = ref_fe->weight_size();
  const size_t dim = ref_fe->dim();

  std::vector< Mat<DIM, DIM, DataType> > fe_weights (fe_weight_size);
  this->hessian_N_fe (ref_pt, fe_ind, fe_weights);

  assert (hessians.size() == dim);
  for (size_t i=0; i<dim; ++i)
  {
    hessians[i] = fe_weights[ref_fe->iv2ind(i,comp)];
  }
}

template<class DataType, int DIM>
void Element<DataType, DIM>::N_fe (const Coord &ref_pt, 
                                   size_t fe_ind, 
                                   std::vector< DataType > &weight) const
{
  doffem::ConstCellTrafoPtr<DataType, DIM> cell_trafo = this->get_cell_transformation();
  assert (cell_trafo->contains_reference_point(ref_pt));
  
  // evaluate shape functions on reference cell
  doffem::RefElement<DataType, DIM> const * ref_fe = this->get_fe(fe_ind);
  const size_t dim = ref_fe->dim();
  const size_t nb_comp = ref_fe->nb_comp();

  assert (weight.size() == ref_fe->weight_size());

  std::vector<DataType> shape_vals (ref_fe->weight_size(), 0.);
  ref_fe->N(ref_pt, shape_vals);

  // map shape function values to element living on physical cell
  IndexFunction ind_fun = boost::bind ( &doffem::RefElement<DataType, DIM>::iv2ind, ref_fe, _1, _2);
  ref_fe->fe_trafo()->map_shape_function_values (*cell_trafo, ref_pt,
                                                 0, dim, nb_comp, ind_fun,
                                                 shape_vals, weight);
}

template<class DataType, int DIM>
void Element<DataType, DIM>::grad_N_fe (const Coord &ref_pt, 
                                        size_t fe_ind, 
                                        std::vector< Vec<DIM, DataType> > &gradients) const
{
  doffem::ConstCellTrafoPtr<DataType, DIM> cell_trafo = this->get_cell_transformation();
  assert (cell_trafo->contains_reference_point(ref_pt));

  // evaluate shape functions values and derivatives on reference cell
  doffem::RefElement<DataType, DIM> const * ref_fe = this->get_fe(fe_ind);
  assert (gradients.size() == ref_fe->weight_size());
  
  const size_t dim = ref_fe->dim();
  const size_t nb_comp = ref_fe->nb_comp();

  std::vector<DataType> shape_vals (ref_fe->weight_size(), 0.);
  std::vector< Vec<DIM,DataType> > shape_grads (ref_fe->weight_size());

  ref_fe->N(ref_pt, shape_vals);
  ref_fe->grad_N(ref_pt, shape_grads);
  
  // evaluate mapped values
  std::vector< DataType > mapped_vals (ref_fe->weight_size());
  this->N_fe (ref_pt, fe_ind, mapped_vals);
  
  // map shape function values to element living on physical cell
  IndexFunction ind_fun = boost::bind ( &doffem::RefElement<DataType, DIM>::iv2ind, ref_fe, _1, _2);
  ref_fe->fe_trafo()->map_shape_function_gradients (*cell_trafo, ref_pt,
                                                    0, dim, nb_comp, ind_fun,
                                                    shape_vals, shape_grads,
                                                    mapped_vals, gradients);
}

template<class DataType, int DIM>
void Element<DataType, DIM>::N_and_grad_N_fe (const Coord &ref_pt, 
                                              size_t fe_ind, 
                                              std::vector< DataType > &weight,
                                              std::vector< Vec<DIM, DataType> > &gradients) const
{
  doffem::ConstCellTrafoPtr<DataType, DIM> cell_trafo = this->get_cell_transformation();
  assert (cell_trafo->contains_reference_point(ref_pt));
  
  // evaluate shape functions on reference cell
  doffem::RefElement<DataType, DIM> const * ref_fe = this->get_fe(fe_ind);
  const size_t dim = ref_fe->dim();
  const size_t nb_comp = ref_fe->nb_comp();

  assert (weight.size() == ref_fe->weight_size());
  assert (gradients.size() == ref_fe->weight_size());

  std::vector<DataType> shape_vals (ref_fe->weight_size(), 0.);
  std::vector< Vec<DIM,DataType> > shape_grads (ref_fe->weight_size());
  
  ref_fe->N(ref_pt, shape_vals);
  ref_fe->grad_N(ref_pt, shape_grads);

  // map shape function values to element living on physical cell
  IndexFunction ind_fun = boost::bind ( &doffem::RefElement<DataType, DIM>::iv2ind, ref_fe, _1, _2);
  ref_fe->fe_trafo()->map_shape_function_values (*cell_trafo, ref_pt,
                                                 0, dim, nb_comp, ind_fun,
                                                 shape_vals, weight);
                                                 

  // map shape function values to element living on physical cell
  ref_fe->fe_trafo()->map_shape_function_gradients (*cell_trafo, ref_pt,
                                                    0, dim, nb_comp, ind_fun,
                                                    shape_vals, shape_grads,
                                                    weight, gradients);
}

template<class DataType, int DIM>
void Element<DataType, DIM>::hessian_N_fe (const Coord &ref_pt, 
                                           size_t fe_ind, 
                                           std::vector< Mat<DIM, DIM, DataType> > &hessians) const
{
  doffem::ConstCellTrafoPtr<DataType, DIM> cell_trafo = this->get_cell_transformation();
  assert (cell_trafo->contains_reference_point(ref_pt));

  // evaluate shape functions values and derivatives on reference cell
  doffem::RefElement<DataType, DIM> const * ref_fe = this->get_fe(fe_ind);
  const size_t dim = ref_fe->dim();
  const size_t nb_comp = ref_fe->nb_comp();

  assert (hessians.size() == ref_fe->weight_size());
  
  std::vector< Vec<DIM,DataType> > shape_grads (ref_fe->weight_size());
  std::vector< Mat<DIM, DIM, DataType> > shape_hessians (ref_fe->weight_size());

  ref_fe->hessian_N(ref_pt, shape_hessians);
  ref_fe->grad_N(ref_pt, shape_grads);

  // evaluate mapped gradients
  std::vector< Vec<DIM, DataType> > mapped_grads (ref_fe->weight_size());
  this->grad_N_fe (ref_pt, fe_ind, mapped_grads);

  // map shape function values to element living on physical cell
  IndexFunction ind_fun = boost::bind ( &doffem::RefElement<DataType, DIM>::iv2ind, ref_fe, _1, _2);
  ref_fe->fe_trafo()->map_shape_function_hessians (*cell_trafo, ref_pt,
                                                   0, dim, nb_comp, ind_fun,
                                                   shape_grads, shape_hessians,   
                                                   mapped_grads, hessians);
}

template class Element<float, 3 >;
template class Element<float, 2 >;
template class Element<float, 1 >;

template class Element<double, 3 >;
template class Element<double, 2 >;
template class Element<double, 1 >;

} // namespace hiflow
