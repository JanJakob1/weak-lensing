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

#ifndef HIFLOW_SPACE_FE_EVALUATION
#define HIFLOW_SPACE_FE_EVALUATION

/// \author Staffan Ronnas, Martin Baumann, Teresa Beck, Philipp Gerstner

#include <map>
#include <set>
#include <string>
#include <vector>

#include "common/vector_algebra.h"
#include "common/parcom.h"
#include "dof/dof_fem_types.h"
#include "mesh/entity.h"
#include <boost/function.hpp>


namespace hiflow {

template <class DataType, int DIM> class VectorSpace;

namespace mesh {
class Entity;
template <class DataType, int DIM> class GeometricSearch;

}

namespace la {
template <class DataType> class Vector;
}

namespace doffem {
template <class DataType, int DIM> class RefElement;
}

template< class DataType, int DIM > 
class FlowTrafo 
{
public:
  FlowTrafo () {}
  
  virtual void operator() (const Vec<DIM, DataType>& coords, 
                           const std::vector<DataType>& flow,
                           std::vector<DataType>& mapped_flow) const = 0;
  
  virtual std::vector<size_t> get_trafo_vars() const = 0;
};

template< class DataType, int DIM > 
class FlowTrafoCyl : public FlowTrafo<DataType, DIM>
{
public:
  FlowTrafoCyl (size_t first_flow_var) 
  {
    trafo_vars_.resize(DIM, 0);
    for (size_t l=0; l<DIM; ++l)
    {
      trafo_vars_[l] = first_flow_var + l;
    }
  }
  
  std::vector<size_t> get_trafo_vars() const
  {
    return this->trafo_vars_;
  }
  
  virtual void operator() (const Vec<DIM, DataType>& coords, 
                           const std::vector<DataType>& flow,
                           std::vector<DataType>& mapped_flow) const
  {
    assert (DIM == 2 || DIM == 3);
    assert (flow.size() == DIM);
    assert (mapped_flow.size() == DIM);
    
    const DataType phi = coords[0];
    mapped_flow[0] = cos(phi) * flow[1] - sin(phi) * flow[0];
    mapped_flow[1] = sin(phi) * flow[1] + cos(phi) * flow[0];
    if (DIM == 3)
    {
      mapped_flow[DIM] = flow[DIM];
    }
  }
  std::vector<size_t> trafo_vars_;
};

template < class DataType, int DIM > 
class FeEvalCell
{
  typedef hiflow::doffem::DofID DofID;
  typedef hiflow::Vec<DIM, DataType> Coord;
  
public:
  /// constructor for evaluating all variables in Fe space for multiple coefficient vectors
  /// return of routine evaluate: vals[iv2ind(0,v)] = v-th comoponent of complete Fe space = v-th variable
  FeEvalCell(const VectorSpace< DataType, DIM > &space, 
             std::vector< hiflow::la::Vector<DataType> const *> coeffs);


  /// constructor for evaluating all variables in Fe space
  /// return of routine evaluate: vals[iv2ind(0,v)] = v-th comoponent of complete Fe space = v-th variable
  FeEvalCell(const VectorSpace< DataType, DIM > &space, 
             const hiflow::la::Vector<DataType> &coeff);
                                 
  /// constructor for evaluating one specific fe 
  /// this type is typically used for FE interpolation
  /// return of routine evaluate: vals[iv2ind(0,v)] = v-th comoponent of specified element
  FeEvalCell(const VectorSpace< DataType, DIM > &space, 
             const hiflow::la::Vector<DataType> &coeff, 
             size_t fe_ind);
                 
  FeEvalCell(const VectorSpace< DataType, DIM > &space, 
             const hiflow::la::Vector<DataType> &coeff, 
             const std::vector<size_t>& fe_ind);
             
  /// constructor for evaluating only specific variables
  /// this type is typically used in cell_visualization
  /// return of routine evaluate: vals[iv2ind(0,v)] = v-th comoponent of complete Fe space,
  /// supposed that v is contained in vars
  FeEvalCell(const VectorSpace< DataType, DIM > &space, 
             const hiflow::la::Vector<DataType> &coeff, 
             const std::vector<size_t>& vars, 
             FlowTrafo<DataType, DIM> const * map);
                 
  virtual ~FeEvalCell() {}

  // TODO: add routines for evaluating at a series of points within the same cell
  void evaluate   (const mesh::Entity& cell, const Coord& pt, std::vector<DataType>& vals) const;
  void r_evaluate (const mesh::Entity& cell, const Coord& ref_pt, std::vector<DataType>& vals) const;
  void r_evaluate (const mesh::Entity& cell, const Coord& ref_pt, std::vector< std::vector<DataType>* > vals) const;
  
  void evaluate_grad   (const mesh::Entity& cell, const Coord& pt, std::vector<Vec<DIM,DataType> >& vals) const;
  void r_evaluate_grad (const mesh::Entity& cell, const Coord& ref_pt, std::vector<Vec<DIM,DataType> >& vals) const;
  void r_evaluate_grad (const mesh::Entity& cell, const Coord& ref_pt, std::vector<std::vector<Vec<DIM,DataType> >* > vals) const;
  
  void r_evaluate_weight_and_grad (const mesh::Entity& cell, 
                                   const Coord& ref_pt, 
                                   std::vector<DataType>& vals,
                                   std::vector<Vec<DIM,DataType> >& grads) const;
                                   
  void r_evaluate_weight_and_grad (const mesh::Entity& cell, 
                                   const Coord& ref_pt, 
                                   std::vector< std::vector<DataType>* > vals,
                                   std::vector< std::vector< Vec<DIM,DataType> >* > grads) const;
                                                            
  inline size_t nb_comp() const
  {
    return this->nb_comp_;
  }

  inline size_t iv2ind(size_t i, size_t v) const
  {
    assert (i==0);
    assert (v >= 0);
    assert (v < this->nb_comp_);
    return v;
  }

  inline size_t ivar2ind(size_t i, size_t v) const
  {
    assert (i==0);
    assert (v >= 0);
    assert (v < this->var_order_.size());
    assert (this->var_order_[v] >= 0);
    return this->var_order_[v];
  }
  
  inline size_t weight_size() const
  {
    return this->nb_comp_;
  }
  
  inline size_t nb_func() const
  {
    return 1;
  }

  void set_print (bool flag)
  {
    this->print_ = flag;
  }
  
protected:
  size_t fe_comp_2_var(size_t fe_ind, size_t comp) const;
 
  void clear_return_values(std::vector<DataType>& vals) const;
  void clear_return_values(std::vector<Vec<DIM,DataType> >& vals) const;

  void update_dof_values(const mesh::Entity& cell) const; 
  
  void setup();
  
  const VectorSpace< DataType, DIM > &space_;
  
  std::vector< la::Vector<DataType> const * > coeffs_;
  
  int gdim_;

  size_t nb_comp_;
  
  std::vector<size_t> vars_;
  std::vector<size_t> fe_inds_;
  std::vector< std::vector< size_t> > fe_ind_2_comp_;
  std::vector< int > var_order_;
  mutable bool print_;
  
  FlowTrafo<DataType, DIM> const * flow_trafo_;
  
  // sorted data sets, i.e. 'val' is sorted corresponding to 'id'
  mutable std::vector< DofID > id_;
  mutable std::vector<std::vector< std::vector<DataType> > > dof_values_; 
  mutable int last_cell_index_;
};


template < class DataType, int DIM > 
class FeEvalLocal
{
  typedef hiflow::doffem::DofID DofID;
  typedef hiflow::Vec<DIM, DataType> Coord;
  
public:
  /// constructor for evaluating all variables in Fe space
  /// return of routine evaluate: vals[iv2ind(0,v)] = v-th comoponent of complete Fe space = v-th variable
  FeEvalLocal(const VectorSpace< DataType, DIM > &space, 
              const hiflow::la::Vector<DataType> &coeff);
                                 
  /// constructor for evaluating one specific fe 
  /// this type is typically used for FE interpolation
  /// return of routine evaluate: vals[iv2ind(0,v)] = v-th comoponent of specified element
  FeEvalLocal(const VectorSpace< DataType, DIM > &space, 
              const hiflow::la::Vector<DataType> &coeff, 
              size_t fe_ind);
                 
  FeEvalLocal(const VectorSpace< DataType, DIM > &space, 
              const hiflow::la::Vector<DataType> &coeff, 
              const std::vector<size_t>& fe_ind);
             
  /// constructor for evaluating only specific variables
  /// this type is typically used in cell_visualization
  /// return of routine evaluate: vals[iv2ind(0,v)] = v-th comoponent of complete Fe space,
  /// supposed that v is contained in vars
  FeEvalLocal(const VectorSpace< DataType, DIM > &space, 
              const hiflow::la::Vector<DataType> &coeff, 
              const std::vector<size_t>& vars, 
              FlowTrafo<DataType, DIM> const * map);
                 
  virtual ~FeEvalLocal();

  void set_trial_cells(const std::vector< int > &trial_cells) const;
  void set_trial_cells(const std::set< int > &trial_cells) const;
    
  bool evaluate (const Coord& pt, DataType& value) const;
  bool evaluate (const Coord& pt, std::vector<DataType>& vals) const;
  
  std::vector<bool> evaluate (const std::vector<Coord>& pt, 
                              std::vector< std::vector<DataType> >& vals) const;

  // here, entity is a dummy argument, for making FeEvalLocal compatible with MappingPhys2Ref
  bool evaluate (const mesh::Entity& entity, const Coord& pt, std::vector<DataType>& vals) const
  {
    return this->evaluate(pt, vals);
  }
  
  // TODO: add routines for evaluating gradient
  bool evaluate_grad              (const Coord& pt, std::vector<Vec<DIM,DataType> > & vals) const;
  std::vector<bool> evaluate_grad (const std::vector<Coord>& pt, 
                                   std::vector< std::vector<Vec<DIM,DataType>> >& vals) const;
       
  inline size_t nb_comp() const
  {
    assert (this->fe_eval_cell_ != nullptr);
    return this->fe_eval_cell_->nb_comp();
  }
  
  inline size_t iv2ind(size_t i, size_t v) const
  {
    assert (this->fe_eval_cell_ != nullptr);
    return this->fe_eval_cell_->iv2ind(i,v);
  }
  
  inline size_t ivar2ind(size_t i, size_t v) const
  {
    assert (this->fe_eval_cell_ != nullptr);
    return this->fe_eval_cell_->ivar2ind(i,v);
  }
  
  inline size_t weight_size() const
  {
    assert (this->fe_eval_cell_ != nullptr);
    return this->fe_eval_cell_->weight_size();
  }
  
  inline size_t nb_func() const
  {
    assert (this->fe_eval_cell_ != nullptr);
    return this->fe_eval_cell_->nb_func();
  }
      
protected:
  virtual std::vector<bool> evaluate_impl (const std::vector<Coord>& pt, 
                                           std::vector< std::vector<DataType> >& vals) const;
                                      
  void setup();
                        
  void search_points (const std::vector<Coord>& pts, 
                      std::vector< std::vector<int> >& cell_indices, 
                      std::vector< std::vector<Coord> >& ref_pts) const;
                      
  bool check_ref_coords (const Coord& pt, 
                         const std::vector<int> & cell_indices, 
                         const std::vector<Coord> & ref_pts) const;
                                               
  mesh::GeometricSearch<DataType, DIM>* search_;
  
  FeEvalCell<DataType, DIM>* fe_eval_cell_;
  
  const VectorSpace< DataType, DIM > &space_;

  mutable bool print_;
  
  mutable std::vector< int > const * vec_trial_cells_;
  
  mutable std::set< int > const * set_trial_cells_;
    
};

template < class DataType, int DIM > 
class FeEvalGlobal : public FeEvalLocal<DataType, DIM>
{
  typedef hiflow::doffem::DofID DofID;
  typedef hiflow::Vec<DIM, DataType> Coord;
  
public:
  /// constructor for evaluating all variables in Fe space
  /// return of routine evaluate: vals[iv2ind(0,v)] = v-th comoponent of complete Fe space = v-th variable
  FeEvalGlobal(const VectorSpace< DataType, DIM > &space, 
               const hiflow::la::Vector<DataType> &coeff);
                                 
  /// constructor for evaluating one specific fe 
  /// this type is typically used for FE interpolation
  /// return of routine evaluate: vals[iv2ind(0,v)] = v-th comoponent of specified element
  FeEvalGlobal(const VectorSpace< DataType, DIM > &space, 
               const hiflow::la::Vector<DataType> &coeff, 
               size_t fe_ind);
                 
  FeEvalGlobal(const VectorSpace< DataType, DIM > &space, 
               const hiflow::la::Vector<DataType> &coeff, 
               const std::vector<size_t>& fe_ind);
             
  /// constructor for evaluating only specific variables
  /// this type is typically used in cell_visualization
  /// return of routine evaluate: vals[iv2ind(0,v)] = v-th comoponent of complete Fe space,
  /// supposed that v is contained in vars
  FeEvalGlobal(const VectorSpace< DataType, DIM > &space, 
               const hiflow::la::Vector<DataType> &coeff, 
               const std::vector<size_t>& vars, 
               FlowTrafo<DataType, DIM> const * map);
                 
  virtual ~FeEvalGlobal();
      
protected:
  std::vector<bool> evaluate_impl (const std::vector<Coord>& pt, 
                                   std::vector< std::vector<DataType> >& vals) const;
                                           
  ParCom* parcom_;
};


template < class DataType, int DIM > 
class BasisEvalLocal
{
  typedef hiflow::Vec<DIM, DataType> Coord;
  typedef hiflow::doffem::DofID BasisId;
  typedef size_t BasisIt;
  typedef size_t CellDofIt;
  typedef int CellIndex;
  typedef size_t CellIt;
  typedef size_t PtIt;
  
public:
  /// constructor for evaluating all basis functions of specified element 
  /// return of routine evaluate: vals[iv2ind(i,v)] = v-th comoponent of the specified fe                                
  BasisEvalLocal(const VectorSpace< DataType, DIM > &space, 
                 size_t fe_ind);

  BasisEvalLocal(const VectorSpace< DataType, DIM > &space, 
                 size_t fe_ind,
                 CellIndex cell_index);

  BasisEvalLocal(const VectorSpace< DataType, DIM > &space, 
                 size_t fe_ind,
                 const std::vector<CellIndex>& cell_indices);

  BasisEvalLocal(const VectorSpace< DataType, DIM > &space, 
                 size_t fe_ind,
                 const std::set<CellIndex>& cell_indices);
                                                                                
  BasisEvalLocal(const VectorSpace< DataType, DIM > &space,
                 size_t fe_ind, bool dummy,
                 const std::vector<BasisId> & global_ids);
                                               
  virtual ~BasisEvalLocal() {}

  void evaluate (const Coord& pt, std::vector<DataType>& vals) const;
  void evaluate (const std::vector<Coord>& pts, 
                 std::vector< std::vector<DataType> >& vals) const;

  // here, entity is a dummy argument, for making BasisEvalLocal compatible with MappingPhys2Ref
  void evaluate (const mesh::Entity& entity, const Coord& pt, std::vector<DataType>& vals) const
  {
    return this->evaluate(pt, vals);
  }
  
  void get_basis_ids (std::vector<BasisId>& basis_ids) const
  {
    basis_ids = this->basis_ids_; 
  }
  
  inline size_t nb_comp() const
  {
    return this->nb_comp_;
  }

  inline size_t iv2ind(size_t i, size_t v) const
  {
    assert (i < this->nb_func_);
    assert (v >= 0);
    assert (v < this->nb_comp_);
    //return i * this->nb_comp_ + v;
    return v * this->nb_func_ + i;
  }

  inline size_t weight_size() const
  {
    return this->weight_size_;
  }
  
  inline size_t nb_func() const
  {
    return this->nb_func_;
  }
    
  void set_print (bool flag)
  {
    this->print_ = flag;
  }
  
protected:
  void clear_return_values(std::vector<DataType>& vals) const;
  void clear_return_values(std::vector<Vec<DIM,DataType> >& vals) const;

  void setup();
  
  void search_pts(const std::vector<Coord>& pts) const; 
  
  const VectorSpace< DataType, DIM > &space_;
   
  int gdim_;

  size_t nb_func_;
  size_t nb_comp_;
  size_t fe_ind_;
  size_t weight_size_;
  mutable bool print_;
  
  mutable const doffem::RefElement< DataType, DIM > * ref_fe_;
  
  // p := point iterator
  // j := basis iterator
  // i := basis id <-> DofPartion global dof id
  // k := cell iterator
  // c := cell index  <-> mesh cell index

  // ...[j] = i =: i[j]
  std::vector<BasisId> basis_ids_;  

  // ...[i] = j, where i = i[j]
  std::map<BasisId, BasisIt> inv_basis_ids_;
    
  // ... [k] = c =: c[k]
  std::vector<CellIndex> active_cells_;
  
  // ... [i] = \{ c: K_{c} \in \supp(phi_{i}) \}
  std::map< BasisId, std::set<CellIndex> > basis_support_;
  
  // ... [c] = \{ i: K_{c} \in \supp(phi_{i}) \}
  std::map< CellIndex, std::set<BasisId> > inv_basis_support_;
  
  // ...[c][i] = dof_factor, with which basis function i, restricted to K_c, has to be multiplied   
  std::map< CellIndex, std::map< BasisId, DataType> > dof_factors_;
  
  // ...[c][i] = l such that \phi_{K_c,l) = \phi_i restricted to K_c_
  std::map< CellIndex, std::map< BasisId, CellDofIt> > global_2_cell_;
  
  // ...[p][i] = #\{ K \in \supp(phi_{i}) : x \in K \}
  //mutable std::vector< std::map< BasisId, int > > pt_multiplicity_;
  
  // ...[p] = #\{ K : x \in K \}
  mutable std::vector< int > pt_multiplicity_;
  
  // ... [p][k] = (pt \in K_{c[k]} ?)
  mutable std::vector< std::vector< bool > > is_pt_in_cell_;
  
  // ... [p] = \{k : p \in K_{c[k]} \}
  mutable std::vector< std::vector< CellIt > > pt_in_cell_;
  
  // ... [p] = \{ref_pt : p \in K_{c[k]}, p <> ref_pt w.r.t. K_{c[k]}\}
  mutable std::vector< std::vector< Coord > > ref_pts_;
  
};

template < class DataType, int DIM > 
class FeEvalBasisLocal
{
  typedef hiflow::Vec<DIM, DataType> Coord;
  typedef hiflow::doffem::DofID DofId;
  typedef size_t BasisIt;
  typedef size_t CellDofIt;
  typedef int CellIndex;
  typedef size_t CellIt;
  typedef size_t PtIt;
  
public:
  /// constructor for evaluating one specific fe 
  /// this type is typically used for FE interpolation
  /// return of routine evaluate: vals[iv2ind(0,v)] = v-th comoponent of specified element
  FeEvalBasisLocal(const VectorSpace< DataType, DIM > &space, 
                   const hiflow::la::Vector<DataType> &coeff, 
                   size_t fe_ind);
                                 
  virtual ~FeEvalBasisLocal();

  //void set_trial_cells(const std::vector< int > &trial_cells) const;
  
  bool evaluate              (const Coord& pt, DataType& value) const;
  bool evaluate              (const Coord& pt, std::vector<DataType>& vals) const;
  std::vector<bool> evaluate (const std::vector<Coord>& pt, 
                              std::vector< std::vector<DataType> >& vals) const;
  
  // here, entity is a dummy argument, for making FeEvalBasisLocal compatible with MappingPhys2Ref
  bool evaluate              (const mesh::Entity& entity, const Coord& pt, std::vector<DataType>& vals) const
  {
    return this->evaluate(pt, vals);
  }
  
  inline size_t nb_comp() const
  {
    return this->nb_comp_;
  }
  
  inline size_t iv2ind(size_t i, size_t v) const
  {
    assert (v < this->nb_comp_);
    assert (i == 0);
    return v;
  }
  
  inline size_t weight_size() const
  {
    return this->weight_size_;
  }
  
  inline size_t nb_func() const
  {
    return this->nb_func_;
  }
      
protected:
  virtual std::vector<bool> evaluate_impl (const std::vector<Coord>& pt, 
                                           std::vector< std::vector<DataType> >& vals) const;
                                      
  BasisEvalLocal<DataType, DIM>* basis_eval_;
  
  const VectorSpace< DataType, DIM > &space_;

  const la::Vector<DataType>& coeff_;
  
  size_t nb_comp_;
  size_t nb_func_;
  size_t weight_size_;
  size_t fe_ind_;
  mutable bool print_;
  mutable std::vector< int > trial_cells_;
};


} // namespace hiflow
#endif
