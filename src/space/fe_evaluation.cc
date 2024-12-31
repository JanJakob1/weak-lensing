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

#include "space/fe_evaluation.h"

#include "common/permutation.h"
#include "fem/fe_manager.h"
#include "fem/fe_reference.h"
#include "fem/cell_trafo/cell_transformation.h"
#include "linear_algebra/vector.h"
#include "mesh/entity.h"
#include "mesh/mesh.h"
#include "mesh/geometric_search.h"
#include "mesh/geometric_tools.h"
#include "space/vector_space.h"
#include "space/element.h"
#include <boost/function.hpp>
#include <set>

namespace hiflow {

template < class DataType >
void sort_dofs (const la::Vector<DataType> &coeff,  
                std::vector< hiflow::doffem::DofID >& id,
                std::vector< DataType >& val)
{
  // get data
  std::vector< hiflow::doffem::DofID > id_tmp;
  std::vector< DataType > val_tmp;

  coeff.GetAllDofsAndValues(id_tmp, val_tmp);

  assert(!id_tmp.empty());
  assert(!val_tmp.empty());
  assert(id_tmp.size() == val_tmp.size());

  // calculate permutation
  std::vector< int > permutation;
  compute_sorting_permutation(id_tmp, permutation);

  // permute
  permute_vector(permutation, id_tmp, id);
  permute_vector(permutation, val_tmp, val);
}

template void sort_dofs <double> (const la::Vector<double> &, std::vector< hiflow::doffem::DofID >&, std::vector< double >&);
template void sort_dofs <float> (const la::Vector<float> &, std::vector< hiflow::doffem::DofID >&, std::vector< float >&);

template <class DataType, int DIM>
void extract_dof_values ( const VectorSpace<DataType, DIM>& space, 
                          const mesh::Entity& cell, 
                          const Element<DataType, DIM>& elem,
                          const la::Vector<DataType>& fun,
                          const std::vector< doffem::DofID >& sorted_id,
                          const std::vector< DataType >& sorted_val,
                          size_t fe_ind, 
                          std::vector<DataType> &dof_values)
{
  typedef hiflow::doffem::DofID DofID;
  std::vector< int > global_dof_ids;
  space.get_dof_indices(fe_ind, cell.index(), &global_dof_ids);

  const int num_dofs = global_dof_ids.size();

  dof_values.clear();
  dof_values.resize(num_dofs, 1.e25);

  size_t num_procs = space.nb_subdom(); 

  if (num_procs == 1) 
  {
    // in sequential world, DofIDs are already sorted
    fun.GetValues(&(global_dof_ids[0]), global_dof_ids.size(), &(dof_values[0]));
  } 
  else 
  {
    // in parallel world, DofIDs are not sorted
    // -> id and val fields need to be sorted and accesses are related to
    //    a seek through the data

    std::vector< DofID >::const_iterator it;
    for (int i = 0; i < num_dofs; ++i) 
    {
      it = std::lower_bound(sorted_id.begin(), sorted_id.end(), global_dof_ids[i]);
      const int index = it - sorted_id.begin();
      dof_values[i] = sorted_val[index];
    }
    // slow version
    // fun_.GetValues(&global_dof_ids[0], num_dofs, &dof_values[0]);
  }
  
  // TODO: is this consistent with previous sorting of dofs?
  std::vector< DataType > dof_factors;
  space.dof().get_dof_factors_on_cell(cell.index(), dof_factors);
  
  const size_t start_dof = elem.dof_offset(fe_ind);
  const size_t end_dof = start_dof + elem.nb_dof(fe_ind);
  for (size_t i=start_dof; i != end_dof; ++i)
  {
    dof_values[i-start_dof] *= dof_factors[i];
  }
}

template void extract_dof_values <float, 1> ( const VectorSpace<float, 1>& , 
                                              const mesh::Entity& ,
                                              const Element<float, 1>& ,
                                              const la::Vector<float>&,
                                              const std::vector< hiflow::doffem::DofID >&,
                                              const std::vector< float >&, 
                                              size_t fe_ind, 
                                              std::vector<float> &);
template void extract_dof_values <float, 2> ( const VectorSpace<float, 2>& , 
                                              const mesh::Entity& ,
                                              const Element<float, 2>& ,
                                              const la::Vector<float>&,
                                              const std::vector< hiflow::doffem::DofID >&,
                                              const std::vector< float >&, 
                                              size_t fe_ind, 
                                              std::vector<float> &);
template void extract_dof_values <float, 3> ( const VectorSpace<float, 3>& , 
                                              const mesh::Entity& ,
                                              const Element<float, 3>& ,
                                              const la::Vector<float>&,
                                              const std::vector< hiflow::doffem::DofID >&,
                                              const std::vector< float >&, 
                                              size_t fe_ind, 
                                              std::vector<float> &);
template void extract_dof_values <double, 1> ( const VectorSpace<double, 1>& , 
                                              const mesh::Entity& , 
                                              const Element<double, 1>& ,
                                              const la::Vector<double>&,
                                              const std::vector< hiflow::doffem::DofID >&,
                                              const std::vector< double >&, 
                                              size_t fe_ind, 
                                              std::vector<double> &);
template void extract_dof_values <double, 2> ( const VectorSpace<double, 2>& , 
                                              const mesh::Entity& , 
                                              const Element<double, 2>& ,
                                              const la::Vector<double>&,
                                              const std::vector< hiflow::doffem::DofID >&,
                                              const std::vector< double >&, 
                                              size_t fe_ind, 
                                              std::vector<double> &);
template void extract_dof_values <double, 3> ( const VectorSpace<double, 3>& , 
                                              const mesh::Entity& , 
                                              const Element<double, 3>& ,
                                              const la::Vector<double>&,
                                              const std::vector< hiflow::doffem::DofID >&,
                                              const std::vector< double >&, 
                                              size_t fe_ind, 
                                              std::vector<double> &);

///////////////////////////////////////////////////////////////////
/////////////// FeEvalCell ////////////////////////////////////
///////////////////////////////////////////////////////////////////
template < class DataType, int DIM >
FeEvalCell<DataType, DIM>::FeEvalCell(const VectorSpace< DataType, DIM > &space, 
                                      std::vector< hiflow::la::Vector<DataType> const *> coeffs)
: space_(space), 
  last_cell_index_(-1),
  flow_trafo_(nullptr)
{
  this->coeffs_ = coeffs;
  this->vars_.clear();
  this->vars_.resize(space.nb_var(), 0);
  for (int l=0; l<this->vars_.size(); ++l)
  {
    this->vars_[l] = l;
  }
  this->setup();
  
}

template < class DataType, int DIM >
FeEvalCell<DataType, DIM>::FeEvalCell(const VectorSpace< DataType, DIM > &space, 
                                      const la::Vector<DataType> &coeff)
: space_(space), 
  last_cell_index_(-1),
  flow_trafo_(nullptr)
{
  this->coeffs_.push_back(&coeff);
  this->vars_.clear();
  this->vars_.resize(space.nb_var(), 0);
  for (int l=0; l<this->vars_.size(); ++l)
  {
    this->vars_[l] = l;
  }
  this->setup();
  
}

template < class DataType, int DIM >
FeEvalCell<DataType, DIM>::FeEvalCell(const VectorSpace< DataType, DIM > &space, 
                                      const la::Vector<DataType> &coeff, 
                                      size_t fe_ind )
: space_(space), 
  last_cell_index_(-1),
  flow_trafo_(nullptr)
{
  assert (fe_ind < space.nb_fe());
  this->coeffs_.push_back(&coeff);
  this->vars_ = space.fe_2_var(fe_ind);
  this->setup();
}

template < class DataType, int DIM >
FeEvalCell<DataType, DIM>::FeEvalCell(const VectorSpace< DataType, DIM > &space, 
                                      const la::Vector<DataType> &coeff, 
                                      const std::vector<size_t> &fe_ind )
: space_(space), 
  last_cell_index_(-1),
  flow_trafo_(nullptr)
{
  this->coeffs_.push_back(&coeff);
  for (size_t l=0; l<fe_ind.size(); ++l)
  {
    assert (fe_ind[l] < space.nb_fe());
    std::vector<size_t> vars = space.fe_2_var(fe_ind[l]);
    
    for (size_t d=0; d<vars.size(); ++d)
    {
      this->vars_.push_back(vars[d]);
    }
  }
  this->setup();
}

template < class DataType, int DIM >
FeEvalCell<DataType, DIM>::FeEvalCell(const VectorSpace< DataType, DIM > &space, 
                                              const la::Vector<DataType> &coeff, 
                                              const std::vector<size_t>& vars,
                                              FlowTrafo<DataType, DIM> const * map )
: space_(space), 
  vars_(vars),
  last_cell_index_(-1),
  flow_trafo_(map)
{
  this->coeffs_.push_back(&coeff);
  this->setup();
}

template < class DataType, int DIM >
void FeEvalCell<DataType, DIM>::setup()
{
  this->gdim_ = this->space_.mesh().gdim();
  const size_t nb_fe = this->space_.nb_fe();
  const size_t nb_var = this->space_.nb_var();
  
  // get fe indices to be evaluated
  this->fe_inds_.clear();
  this->fe_ind_2_comp_.clear();
  this->fe_ind_2_comp_.resize(nb_fe);
  
  std::sort(this->vars_.begin(), this->vars_.end()); 
  auto vlast = std::unique(this->vars_.begin(), this->vars_.end()); 
  this->vars_.erase(vlast, this->vars_.end()); 
  
  this->nb_comp_ = this->vars_.size();
  this->var_order_.clear();
  this->var_order_.resize(nb_var, -1);
  
  for (size_t l=0; l<this->vars_.size(); ++l)
  {
    const size_t v = this->vars_[l];
    assert (v < this->space_.nb_var());
    
    const size_t f = this->space_.var_2_fe(v);
    const size_t c = this->space_.var_2_comp(v);
    
    this->fe_inds_.push_back(f);
    this->fe_ind_2_comp_[f].push_back(c);
    this->var_order_[v] = l;
  }
  
  std::sort(this->fe_inds_.begin(), this->fe_inds_.end()); 
  auto last = std::unique(this->fe_inds_.begin(), this->fe_inds_.end()); 
  this->fe_inds_.erase(last, this->fe_inds_.end()); 
  
  for (size_t f=0; f<nb_fe; ++f)
  {
    if (this->fe_ind_2_comp_[f].size() > 0)
    {
      std::sort(this->fe_ind_2_comp_[f].begin(), this->fe_ind_2_comp_[f].end()); 
      auto flast = std::unique(this->fe_ind_2_comp_[f].begin(), this->fe_ind_2_comp_[f].end()); 
      this->fe_ind_2_comp_[f].erase(flast, this->fe_ind_2_comp_[f].end()); 
    }
  }
  
  const size_t num_coeff = this->coeffs_.size();
  this->dof_values_.clear();
  this->dof_values_.resize(num_coeff);
  for (size_t l=0; l<num_coeff; ++l)
  {
    this->dof_values_[l].resize(nb_fe);
  }
}

template < class DataType, int DIM >
size_t FeEvalCell<DataType, DIM>::fe_comp_2_var(size_t fe_ind, size_t comp) const
{
  return this->space_.fe_2_var(fe_ind)[comp];
}
  
template < class DataType, int DIM >
void FeEvalCell<DataType, DIM>::update_dof_values(const mesh::Entity& cell) const 
{
  int cur_cell_index = cell.index();
  if (cur_cell_index == this->last_cell_index_)
  {
    return;
  }
  
  // TODO: need sort_dofs and extract_dofs as defined above?
  const size_t num_coeff = this->coeffs_.size();
  
  for (size_t i=0; i<num_coeff; ++i)
  {
    for (size_t l=0; l<this->fe_inds_.size(); ++l)
    {
      const size_t f = this->fe_inds_[l];
      //this->space_.set_print(this->print_);
      this->space_.extract_dof_values (f, cur_cell_index, *(this->coeffs_[i]), this->dof_values_[i][f]);
    }
  }
}

template < class DataType, int DIM >
void FeEvalCell<DataType, DIM>::clear_return_values(std::vector<DataType>& vals) const 
{
  if (vals.size() != this->weight_size())
  {
    vals.clear();
    vals.resize(this->weight_size(), 0.);
  }
  else
  {
    for (size_t i=0; i<vals.size(); ++i)
    {
      vals[i] = 0.;
    }
  }
}

template < class DataType, int DIM >
void FeEvalCell<DataType, DIM>::clear_return_values(std::vector<Vec<DIM,DataType> >& vals) const 
{
  if (vals.size() != this->weight_size())
  {
    vals.clear();
    vals.resize(this->weight_size());
  }
  else
  {
    for (size_t i=0; i<vals.size(); ++i)
    {
      vals[i] = Vec<DIM, DataType>();
    }
  }
}

template < class DataType, int DIM >
void FeEvalCell<DataType, DIM>::evaluate (const mesh::Entity& cell, 
                                          const Coord& pt, 
                                          std::vector<DataType>& vals) const
{
  assert (this->coeffs_.size() == 1);
  
  std::shared_ptr<doffem::CellTransformation<DataType, DIM> const > cell_trafo 
    = this->space_.fe_manager().get_cell_transformation(cell.index());
    
  assert (cell_trafo != nullptr);
  
  Coord ref_pt;

  if (!cell_trafo->inverse(pt, ref_pt))
  {
    this->clear_return_values(vals);
  }
  else
  {
    this->r_evaluate (cell, ref_pt, vals);
  }
}

template < class DataType, int DIM >
void FeEvalCell<DataType, DIM>::r_evaluate (const mesh::Entity& cell, 
                                            const Coord& ref_pt, 
                                            std::vector<DataType>& vals) const
{
  assert (this->coeffs_.size() == 1);
  std::vector< std::vector<DataType>* > tmp_vals;
  tmp_vals.push_back(&vals);
  this->r_evaluate(cell, ref_pt, tmp_vals);
}

template < class DataType, int DIM >
void FeEvalCell<DataType, DIM>::r_evaluate_grad (const mesh::Entity& cell, 
                                                 const Coord& ref_pt, 
                                                 std::vector<Vec<DIM,DataType> >& vals) const
{
  assert (this->coeffs_.size() == 1);
  std::vector< std::vector<Vec<DIM,DataType> >* > tmp_vals;
  tmp_vals.push_back(&vals);
  this->r_evaluate_grad(cell, ref_pt, tmp_vals);
}

template < class DataType, int DIM >
void FeEvalCell<DataType, DIM>::r_evaluate_weight_and_grad (const mesh::Entity& cell, 
                                                            const Coord& ref_pt, 
                                                            std::vector<DataType>& vals,
                                                            std::vector<Vec<DIM,DataType> >& grads) const
{
  assert (this->coeffs_.size() == 1);
  
  std::vector< std::vector<DataType>* > tmp_vals;
  tmp_vals.push_back(&vals);
  
  std::vector< std::vector<Vec<DIM,DataType> >* > tmp_grads;
  tmp_grads.push_back(&grads);
  
  this->r_evaluate_weight_and_grad(cell, ref_pt, tmp_vals, tmp_grads);
}

template < class DataType, int DIM >
void FeEvalCell<DataType, DIM>::evaluate_grad (const mesh::Entity& cell, 
                                               const Coord& pt, 
                                               std::vector< Vec<DIM,DataType> >& vals) const
{
  assert (this->coeffs_.size() == 1);

  std::shared_ptr<doffem::CellTransformation<DataType, DIM> const > cell_trafo 
    = this->space_.fe_manager().get_cell_transformation(cell.index());
        
  assert (cell_trafo != nullptr);
  Coord ref_pt;
  if (!cell_trafo->inverse(pt, ref_pt))
  {
    this->clear_return_values(vals);
    return;
  }
  else
  {
    this->r_evaluate_grad (cell, ref_pt, vals);
    return;
  }
}

template < class DataType, int DIM >
void FeEvalCell<DataType, DIM>::r_evaluate (const mesh::Entity& cell, 
                                            const Coord& ref_pt, 
                                            std::vector< std::vector<DataType>* > vals) const
{
  const size_t num_coeff = this->coeffs_.size();
  const size_t num_eval = vals.size();
  assert (num_eval <= num_coeff);

  for (size_t i=0; i<num_eval; ++i)
  {
    assert (vals[i] != nullptr);
    this->clear_return_values(*(vals[i]));
  }
  this->update_dof_values(cell);
     
  Element<DataType, DIM> elem (this->space_, cell.index());

  for (size_t l=0; l<this->fe_inds_.size(); ++l)
  {
    const size_t f = this->fe_inds_[l];
    const size_t dim = elem.get_fe(f)->dim();
    const size_t fe_nb_comp = elem.get_fe(f)->nb_comp();
    const size_t eval_nb_comp = this->fe_ind_2_comp_[f].size();
    
    const size_t weight_size = elem.get_fe(f)->weight_size();
  
    // evaluate mapped FE basis functions in given point
    std::vector<DataType> weights (weight_size, 0.);
    elem.N_fe(ref_pt, f, weights);
    
    for (size_t j=0; j<num_eval; ++j)
    {
      // Global DoF Ids on the given mesh cell
      const size_t num_dofs = this->dof_values_[j][f].size();
      assert (num_dofs == dim);
  
      for (size_t k=0; k<eval_nb_comp; ++k)
      {
        const size_t c = this->fe_ind_2_comp_[f][k];
        const size_t v = this->fe_comp_2_var(f, c);
        const size_t out_index = this->ivar2ind(0,v);
        
        // Summation over weights multiplied by dof_values
        for (size_t i_loc = 0; i_loc < num_dofs; ++i_loc) 
        {      
          vals[j]->at(out_index) += this->dof_values_[j][f][i_loc] * weights[elem.get_fe(f)->iv2ind(i_loc, c)];
        }
      }
    }
  }
  if (this->flow_trafo_ != nullptr)
  {
    const std::vector<size_t> trafo_vars = this->flow_trafo_->get_trafo_vars();
    const size_t nb_trafo_vars = trafo_vars.size();
    std::vector<DataType> trafo_input (nb_trafo_vars);
    std::vector<DataType> trafo_output (nb_trafo_vars);
    
    Vec<DIM, DataType> phys_pt;
    
    std::shared_ptr<doffem::CellTransformation<DataType, DIM> const > cell_trafo 
    = this->space_.fe_manager().get_cell_transformation(cell.index());
        
    cell_trafo->transform(ref_pt, phys_pt);
    
    for (size_t j=0; j<num_eval; ++j)
    {
      for (size_t l=0; l<nb_trafo_vars; ++l)
      { 
        trafo_input[l] = vals[j]->at(this->ivar2ind(0,trafo_vars[l]));
      }
        
      this->flow_trafo_->operator()(phys_pt, trafo_input, trafo_output);
    
      for (size_t l=0; l<nb_trafo_vars; ++l)
      { 
        vals[j]->at(this->ivar2ind(0,trafo_vars[l])) = trafo_output[l];
      }
    }
  }
}
    
template < class DataType, int DIM >
void FeEvalCell<DataType, DIM>::r_evaluate_grad (const mesh::Entity& cell, 
                                                 const Coord& ref_pt, 
                                                 std::vector< std::vector< Vec<DIM,DataType> >* > vals) const
{
  const size_t num_coeff = this->coeffs_.size();
  const size_t num_eval = vals.size();
  assert (num_eval <= num_coeff);

  for (size_t i=0; i<num_eval; ++i)
  {
    assert (vals[i] != nullptr);
    this->clear_return_values(*(vals[i]));
  }
  this->update_dof_values(cell);
     
  Element<DataType, DIM> elem (this->space_, cell.index());

  for (size_t l=0; l<this->fe_inds_.size(); ++l)
  {
    const size_t f = this->fe_inds_[l];
    const size_t dim = elem.get_fe(f)->dim();
    const size_t fe_nb_comp = elem.get_fe(f)->nb_comp();
    const size_t eval_nb_comp = this->fe_ind_2_comp_[f].size();
    
    const size_t weight_size = elem.get_fe(f)->weight_size();
  
    // evaluate mapped FE basis functions in given point
    std::vector<Vec<DIM,DataType> > weights (weight_size);
    elem.grad_N_fe(ref_pt, f, weights);
    
    for (size_t j=0; j<num_eval; ++j)
    {
      // Global DoF Ids on the given mesh cell
      const size_t num_dofs = this->dof_values_[j][f].size();
      assert (num_dofs == dim);
  
      for (size_t k=0; k<eval_nb_comp; ++k)
      {
        const size_t c = this->fe_ind_2_comp_[f][k];
        const size_t v = this->fe_comp_2_var(f, c);
        const size_t out_index = this->ivar2ind(0,v);
        
        // Summation over weights multiplied by dof_values
        for (size_t i_loc = 0; i_loc < num_dofs; ++i_loc) 
        {      
          vals[j]->at(out_index) += this->dof_values_[j][f][i_loc] * weights[elem.get_fe(f)->iv2ind(i_loc, c)];
        }
      }
    }
  }
}

template < class DataType, int DIM >
void FeEvalCell<DataType, DIM>::r_evaluate_weight_and_grad (const mesh::Entity& cell, 
                                                            const Coord& ref_pt, 
                                                            std::vector< std::vector<DataType>* > vals,
                                                            std::vector< std::vector< Vec<DIM,DataType> >* > gradients) const
{
  const size_t num_coeff = this->coeffs_.size();
  const size_t num_eval_fe = vals.size();
  const size_t num_eval_grad = gradients.size();
  
  assert (num_eval_fe <= num_coeff);
  assert (num_eval_grad <= num_coeff);

  assert (num_eval_fe >= num_eval_grad);
  
  for (size_t i=0; i<num_eval_fe; ++i)
  {
    assert (vals[i] != nullptr);
    this->clear_return_values(*(vals[i]));
  }
  
  for (size_t i=0; i<num_eval_grad; ++i)
  {
    assert (gradients[i] != nullptr);
    this->clear_return_values(*(gradients[i]));
  }
  
  this->update_dof_values(cell);
     
  Element<DataType, DIM> elem (this->space_, cell.index());

  for (size_t l=0; l<this->fe_inds_.size(); ++l)
  {
    const size_t f = this->fe_inds_[l];
    const size_t dim = elem.get_fe(f)->dim();
    const size_t fe_nb_comp = elem.get_fe(f)->nb_comp();
    const size_t eval_nb_comp = this->fe_ind_2_comp_[f].size();
    
    const size_t weight_size = elem.get_fe(f)->weight_size();
  
    // evaluate mapped FE basis functions in given point
    std::vector<DataType> basis_vals (weight_size, 0.);
    std::vector<Vec<DIM,DataType> > basis_grads (weight_size);
    
    elem.N_and_grad_N_fe(ref_pt, f, basis_vals, basis_grads);
    
    for (size_t j=0; j<num_eval_fe; ++j)
    {
      // Global DoF Ids on the given mesh cell
      const size_t num_dofs = this->dof_values_[j][f].size();
      assert (num_dofs == dim);
  
      for (size_t k=0; k<eval_nb_comp; ++k)
      {
        const size_t c = this->fe_ind_2_comp_[f][k];
        const size_t v = this->fe_comp_2_var(f, c);
        const size_t out_index = this->ivar2ind(0,v);
        
        // Summation over weights multiplied by dof_values
        for (size_t i_loc = 0; i_loc < num_dofs; ++i_loc) 
        {
          vals[j]->at(out_index) += this->dof_values_[j][f][i_loc] * basis_vals[elem.get_fe(f)->iv2ind(i_loc, c)];      
        }
        
        if (j < num_eval_grad)
        {
          for (size_t i_loc = 0; i_loc < num_dofs; ++i_loc) 
          {
            gradients[j]->at(out_index) += this->dof_values_[j][f][i_loc] * basis_grads[elem.get_fe(f)->iv2ind(i_loc, c)];
          }
        }
      }
    }
  }
}

template class FeEvalCell <float, 1>;
template class FeEvalCell <float, 2>;
template class FeEvalCell <float, 3>;
template class FeEvalCell <double, 1>;
template class FeEvalCell <double, 2>;
template class FeEvalCell <double, 3>;

///////////////////////////////////////////////////////////////////
/////////////// FeEvalLocal ///////////////////////////////////////
///////////////////////////////////////////////////////////////////


template < class DataType, int DIM >
FeEvalLocal<DataType, DIM>::FeEvalLocal(const VectorSpace< DataType, DIM > &space, 
                                        const la::Vector<DataType> &coeff)
: space_(space), search_(nullptr)
{
  this->fe_eval_cell_ = new FeEvalCell<DataType, DIM>(space, coeff);
  this->setup();
}

template < class DataType, int DIM >
FeEvalLocal<DataType, DIM>::FeEvalLocal(const VectorSpace< DataType, DIM > &space, 
                                        const la::Vector<DataType> &coeff, 
                                        size_t fe_ind )
: space_(space), search_(nullptr)
{
  this->fe_eval_cell_ = new FeEvalCell<DataType, DIM>(space, coeff, fe_ind);
  this->setup();
}

template < class DataType, int DIM >
FeEvalLocal<DataType, DIM>::FeEvalLocal(const VectorSpace< DataType, DIM > &space, 
                                        const la::Vector<DataType> &coeff, 
                                        const std::vector<size_t> &fe_ind )
: space_(space), search_(nullptr)
{
  this->fe_eval_cell_ = new FeEvalCell<DataType, DIM>(space, coeff, fe_ind);
  this->setup();
}

template < class DataType, int DIM >
FeEvalLocal<DataType, DIM>::FeEvalLocal(const VectorSpace< DataType, DIM > &space, 
                                        const la::Vector<DataType> &coeff, 
                                        const std::vector<size_t>& vars,
                                        FlowTrafo<DataType, DIM> const * map )
: space_(space), search_(nullptr)
{
  this->fe_eval_cell_ = new FeEvalCell<DataType, DIM>(space, coeff, vars, map);
  this->setup();
}

template < class DataType, int DIM >
FeEvalLocal<DataType, DIM>::~FeEvalLocal()
{
  if (this->fe_eval_cell_ != nullptr)
  {
    delete this->fe_eval_cell_;
  }
  if (this->search_ != nullptr)
  {
    delete this->search_;
  }
}

template < class DataType, int DIM >
void FeEvalLocal< DataType, DIM >::set_trial_cells(const std::vector< int > &trial_cells) const
{
  this->vec_trial_cells_ = &trial_cells;
}

template < class DataType, int DIM >
void FeEvalLocal< DataType, DIM >::set_trial_cells(const std::set< int > &trial_cells) const
{
  this->set_trial_cells_ = &trial_cells;
}

template < class DataType, int DIM >
void FeEvalLocal<DataType, DIM>::setup()
{
  mesh::MeshPtr meshptr = this->space_.meshPtr();
  assert(meshptr != nullptr);

  mesh::GDim gdim = meshptr->gdim();

  if (this->search_ != nullptr)
  {
    delete this->search_;
  }
  
  if (meshptr->is_rectangular()) 
  {
    this->search_ = new mesh::RecGridGeometricSearch<DataType, DIM>(meshptr);
  }   
  else 
  {
    this->search_ = new mesh::GridGeometricSearch<DataType, DIM>(meshptr);
  }
  
  this->vec_trial_cells_ = nullptr;
  this->set_trial_cells_ = nullptr;
}

template < class DataType, int DIM >
void FeEvalLocal<DataType, DIM>::search_points(const std::vector<Coord>& pts, 
                                               std::vector< std::vector<int> >& cell_indices, 
                                               std::vector< std::vector<Coord> >& ref_pts) const
{
  assert (this->search_ != nullptr);
  const size_t nb_pts = pts.size();
  
  if (cell_indices.size() != nb_pts)
  {
    cell_indices.resize(nb_pts);
  }
  if (ref_pts.size() != nb_pts)
  {
    ref_pts.resize(nb_pts);
  }

  for (size_t i=0; i<nb_pts; ++i)
  {
    cell_indices[i].clear();
    ref_pts[i].clear();
  
    if (this->vec_trial_cells_ != nullptr)
    {
      if (!this->vec_trial_cells_->empty()) 
      {
        this->search_->find_cell(pts[i], *this->vec_trial_cells_, cell_indices[i] , ref_pts[i]);
      }
    }
    else if (this->set_trial_cells_ != nullptr)
    {
      if (!this->set_trial_cells_->empty()) 
      {
        this->search_->find_cell(pts[i], *this->set_trial_cells_, cell_indices[i] , ref_pts[i]);
      }
    }
    else 
    {
      this->search_->find_cell(pts[i], cell_indices[i], ref_pts[i]);
    }
#ifndef NDEBUG
    bool success = this->check_ref_coords(pts[i], cell_indices[i], ref_pts[i]);
    assert (success);
#endif
  }
}

template < class DataType, int DIM >
bool FeEvalLocal<DataType, DIM>::evaluate ( const Coord& pt, 
                                            DataType& value ) const 
{
  assert (this->weight_size() == 1);
  std::vector< std::vector< DataType> > tmp_val;
  std::vector< Coord > tmp_pt (1, pt);
  std::vector<bool> found = this->evaluate_impl (tmp_pt, tmp_val);
  assert (tmp_val.size() == 1);
  value = tmp_val[0][0];
  return found[0];
}

template < class DataType, int DIM >
bool FeEvalLocal<DataType, DIM>::evaluate ( const Coord& pt, 
                                            std::vector< DataType >& vals ) const 
{
  std::vector< std::vector< DataType> > tmp_val;
  std::vector< Coord > tmp_pt (1, pt);
  std::vector<bool> found = this->evaluate_impl (tmp_pt, tmp_val);
  assert (tmp_val.size() == 1);
  vals = tmp_val[0];
  return found[0];
}

template < class DataType, int DIM >
std::vector<bool> FeEvalLocal<DataType, DIM>::evaluate ( const std::vector<Coord>& pts, 
                                                         std::vector<std::vector<DataType> >& vals ) const 
{
  return this->evaluate_impl(pts, vals); 
}

template < class DataType, int DIM >
std::vector<bool> FeEvalLocal<DataType, DIM>::evaluate_impl ( const std::vector<Coord>& pts, 
                                                              std::vector<std::vector<DataType> >& vals ) const 
{
  assert (this->search_ != nullptr);
  std::vector<bool> success (pts.size(), true);
  
  if (vals.size() != pts.size())
  {
    vals.resize(pts.size());
  }
  
  mesh::MeshPtr meshptr = this->space_.meshPtr();
  
  const size_t w_size = this->weight_size();
  
  // get cell indices and reference coords for given physical points
  std::vector< std::vector<int> > cell_indices; 
  std::vector< std::vector<Coord> > ref_pts;
  this->search_points(pts, cell_indices, ref_pts);
  
  assert (ref_pts.size() == cell_indices.size());
  
  for (size_t p=0; p<pts.size(); ++p)
  {
    int value_count = ref_pts[p].size();
    vals[p].clear();
    vals[p].resize(w_size, 0.);
    assert (ref_pts[p].size() == cell_indices[p].size());
    
    /*
    if (print_)
    {
      std::cout << "FeEvalLocal: phys point " << pts[p] << " multipl. " << value_count << std::endl;
    }
    */
    
    if (value_count > 0) 
    {
      // point was found in local cells
      for (size_t i = 0; i < value_count; ++i) 
      {
        mesh::Entity cell = meshptr->get_entity(meshptr->tdim(), cell_indices[p][i]);
        std::vector<DataType> tmp_vals (w_size, 0.);
        
        this->fe_eval_cell_->set_print(this->print_);
        this->fe_eval_cell_->r_evaluate (cell, ref_pts[p][i], tmp_vals);
        
        for (size_t l=0; l<w_size; ++l)
        {
          vals[p][l] += tmp_vals[l];
        }
      } 
      for (size_t l=0; l<w_size; ++l)
      {
        vals[p][l] *= 1. / static_cast< DataType >(value_count);
      }
      
      //std::cout << "     ->  " << string_from_range(vals[p].begin(), vals[p].end()) << std::endl;
    }
    else
    {
      success[p] = false;
    }
  }
  return success;
}


template < class DataType, int DIM >
bool FeEvalLocal<DataType, DIM>::check_ref_coords(const Coord& pt, 
                                                  const std::vector<int> & cell_indices, 
                                                  const std::vector<Coord> & ref_pts) const
{
  const DataType eps = 1e-8;
  assert (cell_indices.size() == ref_pts.size());
  
  for (size_t i=0; i<cell_indices.size(); ++i)
  {
    std::shared_ptr<const doffem::CellTransformation<DataType, DIM>> trafo =
       this->space_.get_cell_transformation(cell_indices[i]);
  
    Coord mpt;
    trafo->transform(ref_pts[i], mpt);
    
    DataType diff = norm(mpt-pt);
    if (diff > eps)
    {
      // std::cout << trafo->name() << std::endl;
      //std::cout << pt << " <=> " << ref_pts[i] << " <=> " << mpt << std::endl;
      //trafo->print_vertex_coords();
      return false;
    }
  }
  return true;
}

template class FeEvalLocal <float, 1>;
template class FeEvalLocal <float, 2>;
template class FeEvalLocal <float, 3>;
template class FeEvalLocal <double, 1>;
template class FeEvalLocal <double, 2>;
template class FeEvalLocal <double, 3>;

///////////////////////////////////////////////////////////////////
/////////////// FeEvalGlobal //////////////////////////////////////
///////////////////////////////////////////////////////////////////


template < class DataType, int DIM >
FeEvalGlobal<DataType, DIM>::FeEvalGlobal(const VectorSpace< DataType, DIM > &space, 
                                          const la::Vector<DataType> &coeff)
: FeEvalLocal<DataType, DIM> (space, coeff)
{
  this->parcom_ = new ParCom(this->space_.get_mpi_comm());
}

template < class DataType, int DIM >
FeEvalGlobal<DataType, DIM>::FeEvalGlobal(const VectorSpace< DataType, DIM > &space, 
                                          const la::Vector<DataType> &coeff, 
                                          size_t fe_ind )
: FeEvalLocal<DataType, DIM> (space, coeff, fe_ind)
{
  this->parcom_ = new ParCom(this->space_.get_mpi_comm());
}

template < class DataType, int DIM >
FeEvalGlobal<DataType, DIM>::FeEvalGlobal(const VectorSpace< DataType, DIM > &space, 
                                          const la::Vector<DataType> &coeff, 
                                          const std::vector<size_t> &fe_ind )
: FeEvalLocal<DataType, DIM> (space, coeff, fe_ind)
{
  this->parcom_ = new ParCom(this->space_.get_mpi_comm());
}

template < class DataType, int DIM >
FeEvalGlobal<DataType, DIM>::FeEvalGlobal(const VectorSpace< DataType, DIM > &space, 
                                          const la::Vector<DataType> &coeff, 
                                          const std::vector<size_t>& vars,
                                          FlowTrafo<DataType, DIM> const * map )
: FeEvalLocal<DataType, DIM> (space, coeff, vars, map)
{
  this->parcom_ = new ParCom(this->space_.get_mpi_comm());
}

template < class DataType, int DIM >
FeEvalGlobal<DataType, DIM>::~FeEvalGlobal()
{
  if (this->parcom_ != nullptr)
  {
    delete this->parcom_;
  }
}

template < class DataType, int DIM >
std::vector<bool> FeEvalGlobal<DataType, DIM>::evaluate_impl (const std::vector<Coord>& pts, 
                                                              std::vector<std::vector<DataType> >& vals ) const 
{
  const size_t num_pt = pts.size();
  
  std::vector<bool> success (num_pt, true);
  
  if (vals.size() != num_pt)
  {
    vals.resize(num_pt);
  }

  std::vector<std::vector<DataType> > local_vals;
  std::vector<bool> local_found = FeEvalLocal<DataType, DIM>::evaluate_impl( pts, local_vals );

  assert (local_vals.size() == local_found.size());
  assert (local_found.size() == num_pt);
   
  this->parcom_->sum(local_vals, vals);
   
  assert (local_vals.size() == vals.size());
  
  std::vector<DataType> local_denom (num_pt, 0.);
  std::vector<DataType> denom (num_pt, 0.);
  
  for (size_t d=0; d<num_pt; ++d)
  {
    local_denom[d] = static_cast<DataType>(local_found[d]);
  } 
  
  this->parcom_->sum(local_denom, denom);  

  for (size_t d=0; d<num_pt; ++d)
  {
    if (denom[d] > 0)
    {
      const size_t size_d = vals[d].size();
      for (size_t i=0; i<size_d; ++i)
      {
        vals[d][i] /= denom[d];
      }
    }
    else
    {
      success[d] = false;
    }
  }
  return success;
}

template class FeEvalGlobal <float, 1>;
template class FeEvalGlobal <float, 2>;
template class FeEvalGlobal <float, 3>;
template class FeEvalGlobal <double, 1>;
template class FeEvalGlobal <double, 2>;
template class FeEvalGlobal <double, 3>;

///////////////////////////////////////////////////////////////////
/////////////// BasisEvalLocal ////////////////////////////////////
///////////////////////////////////////////////////////////////////


template < class DataType, int DIM >
BasisEvalLocal<DataType, DIM>::BasisEvalLocal(const VectorSpace< DataType, DIM > &space, 
                                              size_t fe_ind)
: space_(space), fe_ind_(fe_ind)
{
  std::set<BasisId> tmp_ids;
  const int tdim = this->space_.meshPtr()->tdim();
  
  for (mesh::EntityIterator cell = this->space_.meshPtr()->begin(tdim); 
       cell != this->space_.meshPtr()->end(tdim); ++cell) 
  {
    std::vector< BasisId > gl_dofs_on_cell;
    this->space_.get_dof_indices(fe_ind, cell->index(), &gl_dofs_on_cell);
    tmp_ids.insert(gl_dofs_on_cell.begin(), gl_dofs_on_cell.end());
  }
  this->basis_ids_.clear();
  
  auto e_it = tmp_ids.end();
  for (auto it = tmp_ids.begin(); it != e_it; ++it)
  {
    this->basis_ids_.push_back(*it);
  }
  this->setup();
}

template < class DataType, int DIM >
BasisEvalLocal<DataType, DIM>::BasisEvalLocal(const VectorSpace< DataType, DIM > &space, 
                                              size_t fe_ind,
                                              CellIndex cell_index)
: space_(space), fe_ind_(fe_ind)
{
  this->space_.get_dof_indices(this->fe_ind_, cell_index, &this->basis_ids_);
  this->setup();
}

template < class DataType, int DIM >
BasisEvalLocal<DataType, DIM>::BasisEvalLocal(const VectorSpace< DataType, DIM > &space, 
                                              size_t fe_ind,
                                              const std::vector<CellIndex>& cell_indices)
: space_(space), fe_ind_(fe_ind)
{
  std::set<BasisId> all_ids;
  for (size_t l=0; l<cell_indices.size(); ++l)
  {
    const CellIndex c = cell_indices[l];
    std::vector<BasisId> tmp_ids;

    assert (c >= 0);
    assert (c < this->space_.meshPtr()->num_entities(DIM));
    assert (c < this->space_.fe_manager().fe_tank_size());
    
    this->space_.get_dof_indices(this->fe_ind_, c, &tmp_ids);
    all_ids.insert(tmp_ids.begin(), tmp_ids.end());
  }
  
  auto e_it = all_ids.end();
  for (auto it = all_ids.begin(); it != e_it; ++it)
  {
    this->basis_ids_.push_back(*it);
  }
  this->setup();
}

template < class DataType, int DIM >
BasisEvalLocal<DataType, DIM>::BasisEvalLocal(const VectorSpace< DataType, DIM > &space, 
                                              size_t fe_ind,
                                              const std::set<CellIndex>& cell_indices)
: space_(space), fe_ind_(fe_ind)
{
  std::set<BasisId> all_ids;
  auto e_it = cell_indices.end();
  for (auto it = cell_indices.begin(); it != e_it; ++it)
  {
    assert (*it >= 0);
    assert (*it < this->space_.meshPtr()->num_entities(DIM));
    assert (*it < this->space_.fe_manager().fe_tank_size());
    
    std::vector<BasisId> tmp_ids;
    this->space_.get_dof_indices(this->fe_ind_, *it, &tmp_ids);
    all_ids.insert(tmp_ids.begin(), tmp_ids.end());
  }
  
  auto e_it2 = all_ids.end();
  for (auto it = all_ids.begin(); it != e_it2; ++it)
  {
    this->basis_ids_.push_back(*it);
  }
  this->setup();
}

template < class DataType, int DIM >
BasisEvalLocal<DataType, DIM>::BasisEvalLocal(const VectorSpace< DataType, DIM > &space, 
                                              size_t fe_ind, bool dummy,
                                              const std::vector<BasisId>& global_ids )
: space_(space), basis_ids_(global_ids), fe_ind_(fe_ind)
{
  this->setup();
}

template < class DataType, int DIM >
void BasisEvalLocal<DataType, DIM>::setup( )
{
  std::vector<size_t> vars = this->space_.fe_2_var (this->fe_ind_);
  this->nb_func_ = this->basis_ids_.size();
  this->nb_comp_ = vars.size();
  this->weight_size_ = this->nb_func_ * this->nb_comp_;
 
  std::set<CellIndex> cells;
  
  this->inv_basis_ids_.clear();
  this->basis_support_.clear();
  this->inv_basis_support_.clear();
  
  for (BasisIt j=0; j<this->nb_func_; ++j)
  {
    const BasisId i = this->basis_ids_[j];
    this->inv_basis_ids_[i] = j;
    
    std::vector<CellIndex> cur_cells;
    this->space_.dof().global2cells (i, cur_cells);
        
    cells.insert(cur_cells.begin(), cur_cells.end());

    assert (this->basis_support_[i].size() == 0);
    this->basis_support_[i].insert(cur_cells.begin(), cur_cells.end());
    for (size_t l=0; l<cur_cells.size(); ++l)
    {
      this->inv_basis_support_[cur_cells[l]].insert(i);
    }
  }
  
  this->active_cells_.clear();
  auto e_it = cells.end();
  for (auto c_it = cells.begin(); c_it != e_it; ++c_it)
  {
    this->active_cells_.push_back(*c_it);
  }
  
  this->dof_factors_.clear();
  this->global_2_cell_.clear();
      
  for (CellIt k=0; k<this->active_cells_.size(); ++k)
  {
    const CellIndex c = this->active_cells_[k];
    
    // get dof factors
    size_t start_dof = 0;
    for (size_t f=0; f<this->fe_ind_; ++f)
    { 
      start_dof += this->space_.fe_manager().get_fe(c, f)->dim();
    }
  
    std::vector< BasisId > all_gl_dof_ids_on_cell;
    std::vector< DataType> cur_dof_factors;
    
    this->space_.get_dof_indices(this->fe_ind_, c, &all_gl_dof_ids_on_cell);
    this->space_.dof().get_dof_factors_on_cell(c, cur_dof_factors);
    
    const size_t nb_dof_fe = all_gl_dof_ids_on_cell.size();
    for (CellDofIt l = 0; l < nb_dof_fe; ++l)
    {
      this->dof_factors_[c][all_gl_dof_ids_on_cell[l]] = cur_dof_factors[start_dof + l];
      this->global_2_cell_[c][all_gl_dof_ids_on_cell[l]] = l;
    }
  }
}

template < class DataType, int DIM >
void BasisEvalLocal<DataType, DIM>::evaluate (const Coord& pt, std::vector<DataType>& vals) const
{
  std::vector<Coord> tmp_pts(1, pt);
  std::vector< std::vector< DataType> > tmp_vals;
  this->evaluate(tmp_pts, tmp_vals);
  assert (tmp_vals.size() == 1);
  vals = tmp_vals[0];
}

template < class DataType, int DIM >
void BasisEvalLocal<DataType, DIM>::evaluate (const std::vector<Coord>& pts, 
                                              std::vector< std::vector<DataType> >& vals) const
{
  const size_t num_pt = pts.size();
  const size_t num_cell = this->active_cells_.size();
  const int tdim = this->space_.meshPtr()->tdim();
  
  vals.resize(num_pt);
  
  // setup data structures for current points
  this->search_pts(pts);
  
  std::vector<DataType> cur_weights;
  
  // loop over all points
  for (PtIt p=0; p<num_pt; ++p)
  {
    vals[p].clear();
    vals[p].resize(this->weight_size_, 0.);
    //const Coord pt = pts[p];
    const DataType mult_factor = 1. / static_cast<DataType> (this->pt_multiplicity_[p]);
    
    // loop over all active cells (= Union of support of all considered basis functions)
    // that contain the current point
    for (size_t l=0; l<pt_in_cell_[p].size(); ++l)
    {
      // cell K_c
      const CellIt k = this->pt_in_cell_[p][l];
      const CellIndex c = this->active_cells_[k]; 
      const Coord ref_pt = this->ref_pts_[p][l];

      /*
      if (this->print_)
      {
        std::cout << "BasisEval at " << pt << " = ref pt " << ref_pt << " on cell " << c << std::endl; 
      }
      */
      
      assert (this->pt_in_cell_[p].size() == this->ref_pts_[p].size());
      assert (this->is_pt_in_cell_[p][k]);
      assert (this->active_cells_[k] < this->space_.meshPtr()->num_entities(tdim));
            
      // evaluate mapped FE basis on current cell at current point
      Element<DataType, DIM> elem (this->space_, c);
      this->ref_fe_ = elem.get_fe(this->fe_ind_);
          
      const size_t cur_weight_size = elem.get_fe(this->fe_ind_)->weight_size();
      const size_t nb_dof_fe = elem.nb_dof(this->fe_ind_);
      
      cur_weights.clear();
      cur_weights.resize (cur_weight_size, 0.);
      elem.N_fe(ref_pt, this->fe_ind_, cur_weights);
  
      // loop over all basis functions, whose respective support has a non-empty intersection with K_c
      auto end_i_it = this->inv_basis_support_.at(c).end();
      assert (end_i_it != this->inv_basis_support_.at(c).begin());
      
      for (auto i_it = this->inv_basis_support_.at(c).begin(); i_it != end_i_it; ++i_it)
      {
        //assert (this->dof_factors_[c].find(*i_it) != this->dof_factors_[c].end());
        assert (this->global_2_cell_.at(c).find(*i_it) != this->global_2_cell_.at(c).end());
        //assert (this->pt_multiplicity_[p].find(*i_it) != this->pt_multiplicity_[p].end());
        assert (this->inv_basis_ids_.find(*i_it) != this->inv_basis_ids_.end());
        
        const CellDofIt loc_i = this->global_2_cell_.at(c).at(*i_it);
        const BasisIt j = this->inv_basis_ids_.at(*i_it);
        assert (loc_i >= 0);
        assert (loc_i < nb_dof_fe);
        assert (j >= 0);
        assert (j < this->nb_func_);
        const DataType dof_factor = this->dof_factors_.at(c).at(*i_it);
//      const DataType mult_factor = 1. / static_cast<DataType> (this->pt_multiplicity_[p].at(*i_it));

        /*        
        if (this->print_)
        {
          std::cout << " cell " << c << " dof id " << *i_it 
          *         << " dof factor " << dof_factor << " mult_factor " << mult_factor << std::endl;
        }
        */
        
        // loop over components of FE
        for (size_t cp = 0; cp < this->nb_comp_; ++cp)
        {
          vals[p][this->iv2ind(j, cp)] += dof_factor 
                                        * mult_factor 
                                        * cur_weights[this->ref_fe_->iv2ind(loc_i, cp)];
        }
      }
    }
  }
}
  
template < class DataType, int DIM >
void BasisEvalLocal<DataType, DIM>::search_pts(const std::vector<Coord>& pts) const
{
  const size_t num_pt = pts.size();
  const size_t num_cell = this->active_cells_.size();
  const int tdim = this->space_.meshPtr()->tdim();
  
  // todo: for performance reasons: try to avoid clear / resize
  this->pt_multiplicity_.clear();
  this->is_pt_in_cell_.clear();
  this->pt_in_cell_.clear();
  this->ref_pts_.clear();
  
  this->pt_multiplicity_.resize(num_pt, 0);
  this->is_pt_in_cell_.resize(num_pt);
  this->pt_in_cell_.resize(num_pt);
  this->ref_pts_.resize(num_pt);
  
  for (PtIt p = 0; p<num_pt; ++p)
  {
    this->is_pt_in_cell_[p].resize(num_cell, false);
  }
      
  std::vector<DataType> cell_coord;
  std::vector<double> double_coord;
  
  // loop over active cells
  for (CellIt k=0; k<num_cell; ++k)
  {
    const CellIndex c = this->active_cells_[k];
    
    assert (this->inv_basis_support_.find(c) != this->inv_basis_support_.end());
    
    // get vertex coordinates
    double_coord.clear();
    
    this->space_.meshPtr()->get_coordinates(tdim, c, double_coord);
    
    cell_coord.clear();
    
    double_2_datatype(double_coord, cell_coord);
       
    // loop over points
    for (PtIt p = 0; p<num_pt; ++p)
    {
      const Coord pt = pts[p];
      Coord ref_pt;
      
      // check if pt is contained in current cell 
      bool found = mesh::point_inside_cell<DataType, DIM>(pt, cell_coord, ref_pt);

      this->is_pt_in_cell_[p][k] = found;
      if (found)
      {
        this->pt_in_cell_[p].push_back(k);
        this->ref_pts_[p].push_back(ref_pt);
        this->pt_multiplicity_[p]++;
        
/*
        auto e_it = this->inv_basis_support_.at(c).end();
        for (auto i=this->inv_basis_support_.at(c).begin(); i!=e_it; ++i)
        {
          assert (*i >= 0);
          if (this->pt_multiplicity_[p].find(*i) == this->pt_multiplicity_[p].end())
          {
            this->pt_multiplicity_[p][*i] = 1;
          }
          else
          {
            this->pt_multiplicity_[p][*i] += 1;
          }
        }
*/
      }
    }
  }
}

template class BasisEvalLocal <float, 1>;
template class BasisEvalLocal <float, 2>;
template class BasisEvalLocal <float, 3>;
template class BasisEvalLocal <double, 1>;
template class BasisEvalLocal <double, 2>;
template class BasisEvalLocal <double, 3>;


///////////////////////////////////////////////////////////////////
/////////////// FeEvalLocal ///////////////////////////////////////
///////////////////////////////////////////////////////////////////

template < class DataType, int DIM >
FeEvalBasisLocal<DataType, DIM>::FeEvalBasisLocal(const VectorSpace< DataType, DIM > &space, 
                                                  const la::Vector<DataType> &coeff, 
                                                  size_t fe_ind )
: space_(space), fe_ind_(fe_ind), coeff_(coeff)
{
  this->nb_comp_ = space.fe_2_var(fe_ind).size();
  this->nb_func_ = 1;
  this->weight_size_ = this->nb_comp_;
  this->basis_eval_ = new BasisEvalLocal<DataType, DIM>(space, fe_ind);

  assert (this->basis_eval_->nb_comp() == this->nb_comp());
}

template < class DataType, int DIM >
FeEvalBasisLocal<DataType, DIM>::~FeEvalBasisLocal()
{
  if (this->basis_eval_ != nullptr)
  {
    delete this->basis_eval_;
  }
}

template < class DataType, int DIM >
bool FeEvalBasisLocal<DataType, DIM>::evaluate ( const Coord& pt, 
                                                 DataType& value ) const 
{
  assert (this->weight_size() == 1);
  std::vector< std::vector< DataType> > tmp_val;
  std::vector< Coord > tmp_pt (1, pt);
  std::vector<bool> found = this->evaluate_impl (tmp_pt, tmp_val);
  assert (tmp_val.size() == 1);
  value = tmp_val[0][0];
  return found[0];
}

template < class DataType, int DIM >
bool FeEvalBasisLocal<DataType, DIM>::evaluate ( const Coord& pt, 
                                                 std::vector< DataType >& vals ) const 
{
  std::vector< std::vector< DataType> > tmp_val;
  std::vector< Coord > tmp_pt (1, pt);
  std::vector<bool> found = this->evaluate_impl (tmp_pt, tmp_val);
  assert (tmp_val.size() == 1);
  vals = tmp_val[0];
  return found[0];
}

template < class DataType, int DIM >
std::vector<bool> FeEvalBasisLocal<DataType, DIM>::evaluate ( const std::vector<Coord>& pts, 
                                                              std::vector<std::vector<DataType> >& vals ) const 
{
  return this->evaluate_impl(pts, vals); 
}

template < class DataType, int DIM >
std::vector<bool> FeEvalBasisLocal<DataType, DIM>::evaluate_impl ( const std::vector<Coord>& pts, 
                                                                   std::vector<std::vector<DataType> >& vals ) const 
{
  std::vector<bool> success (pts.size(), true);
  
  if (vals.size() != pts.size())
  {
    vals.resize(pts.size());
  }

  /*
  if (this->print_)
  {  
    for (size_t p=0; p<pts.size(); ++p)
    {
      std::cout << "FeEvalBasisLocal at " << pts[p] << std::endl;
    }  
  }
  */
  
  mesh::MeshPtr meshptr = this->space_.meshPtr();
  
  const size_t w_size = this->weight_size();
  
  std::vector< std::vector<DataType> > basis_vals;
  std::vector< DofId > basis_ids;

  //this->basis_eval_->set_print(this->print_);
  this->basis_eval_->evaluate(pts, basis_vals);
  this->basis_eval_->get_basis_ids(basis_ids);
  
  //std::cout << string_from_range(basis_ids.begin(), basis_ids.end()) << std::endl;
  
  const size_t num_basis = basis_ids.size();
  
  std::vector< DataType > coeff_val (num_basis, 0.);
  this->coeff_.GetValues(&basis_ids[0], num_basis, &coeff_val[0]);
  
  for (size_t p=0; p<pts.size(); ++p)
  {
    vals[p].clear();
    vals[p].resize(w_size, 0.);
    
    for (size_t j=0; j<num_basis; ++j)
    {
      for (size_t cp = 0; cp<this->nb_comp_; ++cp)
      {
        vals[p][this->iv2ind(0,cp)] += coeff_val[j] * basis_vals[p][basis_eval_->iv2ind(j,cp)];
      }
    }
  }
  return success;
}


template class FeEvalBasisLocal <float, 1>;
template class FeEvalBasisLocal <float, 2>;
template class FeEvalBasisLocal <float, 3>;
template class FeEvalBasisLocal <double, 1>;
template class FeEvalBasisLocal <double, 2>;
template class FeEvalBasisLocal <double, 3>;


} // namespace hiflow
