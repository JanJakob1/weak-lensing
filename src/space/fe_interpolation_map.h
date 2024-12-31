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

#ifndef HIFLOW_SPACE_FE_INTERPOLATION_MAP
#define HIFLOW_SPACE_FE_INTERPOLATION_MAP

/// \author Philipp Gerstner

#include <map>
#include <vector>

#include "common/data_tools.h"
#include "common/parcom.h"
#include "common/vector_algebra.h"
#include "common/sorted_array.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/lmp/init_vec_mat.h"
#include "linear_algebra/lmp/lmatrix.h"
#include "linear_algebra/lmp/lmatrix_csr_cpu.h"
#include "fem/fe_mapping.h"
#include "mesh/geometric_tools.h"
#include "mesh/entity.h"
#include "space/vector_space.h"
#include "space/fe_evaluation.h"
#include "space/fe_interpolation_cell.h"
#include "space/fe_evaluation.h"

namespace hiflow {

template < class LAD, int DIM > 
class FeInterMapBase
{
  typedef typename LAD::VectorType Vector;
  typedef typename LAD::DataType DataType;
  typedef hiflow::Vec<DIM, DataType> Coord;
  typedef hiflow::doffem::DofID DofId;
  typedef size_t BasisIt;
  typedef size_t CellDofIt;
  typedef int CellIndex;
  typedef size_t CellIt;
  typedef size_t PtIt;
  
public:

  FeInterMapBase();

  ~FeInterMapBase()
  {
    this->clear(); 
  }
 
  inline bool is_initialized() const
  {
    return this->initialized_;
  }
  
protected:
  void clear();

  template <class CellInterpolator>
  void init_with_linear_map_without_comm (CellInterpolator * cell_inter,
                                          VectorSpace< DataType, DIM> const * in_space,
                                          VectorSpace< DataType, DIM> const * out_space,
                                          const Vector& in_vec, 
                                          const Vector& out_vec,
                                          const std::vector<size_t>& in_fe_inds, 
                                          const std::vector<size_t>& out_fe_inds);
                     
  void init_without_linear_map_without_comm (VectorSpace< DataType, DIM> const * in_space,
                                             VectorSpace< DataType, DIM> const * out_space,
                                             const std::vector<size_t>& in_fe_inds, 
                                             const std::vector<size_t>& out_fe_inds);
  
  void build_interpolation_matrix(const Vector& in_vec, const Vector& out_vec);
    
  void interpolate_with_linear_map_without_comm (const Vector& in_vec, 
                                                 Vector& out_vec) const;
                                                     
  template <class CellInterpolator>
  void interpolate_without_linear_map_without_comm (CellInterpolator * cell_inter, 
                                                    const Vector& in_vec, 
                                                    Vector& out_vec) const;
  
  void interpolate_with_matrix_without_comm (const Vector& in_vec, 
                                             Vector& out_vec) const;

  void interpolate_trans_with_matrix_without_comm (const Vector& in_vec, 
                                                   Vector& out_vec) const;
                                                                                                                        
  VectorSpace< DataType, DIM> const * in_space_;
  VectorSpace< DataType, DIM> const * out_space_;
  
  la::lMatrix< DataType > * diag_;
  la::lMatrix< DataType > * diagT_;
  la::lMatrix< DataType > * odiag_;
  la::lMatrix< DataType > * odiagT_;

  //std::vector<DofId> in_ids_;
  //std::vector<DofId> out_ids_;
  
  /// \brief cell_map[K] = {K' \in in_mesh : K \cap K' != \emptyset}
  std::map< int, std::set<int> > cell_map_;
  
  /// \brief weights[K] = { A^e_{K,K'} \in R^{m, n} : K' \in cell_map[K]}
  /// with A^e_{K,K'} [i,j] = [cell_inter.eval(phi_{K', j, e})]_i 
  /// m : #dofs on cell K
  /// n : #dofs on cell K'
  
  //std::vector< std::map< int, std::vector< std::vector< std::vector< DataType > > > > > weights_;
  std::vector< std::vector< std::vector< std::vector< DataType > > > > weights_;
  
  /// \brief in_dof_ids[K][e] = global dof ids of basis functions of FE index e of in_space 
  /// whose support has nonempty intersection with out_cell K 
  //std::vector< std::vector< std::vector< doffem::DofID > > > in_dof_ids_;
  std::vector< std::vector< std::vector< DofId > > > in_dof_ids_;
  
  //std::map < doffem::DofID, DataType > in_dof_weights_;
  
  /// \brief out_dof_ids[K][e] = global dof ids on cell K w.r.t. to FE e in out_space
  std::vector< std::vector< std::vector< DofId > > > out_dof_ids_;
  
  std::vector<size_t> in_fe_inds_;
  std::vector<size_t> out_fe_inds_;

  mutable std::vector< DataType > in_dof_factors_;
  mutable std::vector< DataType > out_dof_factors_;
  mutable std::vector< DataType > in_vals_;
  mutable std::vector< DataType > out_dof_values_;
  
  bool initialized_;
  
  ParCom* parcom_;
};

/// \brief Interpolation map for nodal interpolation
/// For convenience regarding the ugly template argument "CellInterpolator" 
template < class LAD, int DIM> 
class FeInterMapFullNodal : public FeInterMapBase<LAD, DIM> 
{
  typedef typename LAD::DataType DataType;
  typedef typename LAD::VectorType Vector;
  typedef FeInterCellNodal<DataType, DIM, BasisEvalLocal< DataType, DIM > >  NodalCellInterpolator;
  
public:

  FeInterMapFullNodal()
  : FeInterMapBase<LAD, DIM >()
  {
  }

  ~FeInterMapFullNodal()
  {
  }

  void init (VectorSpace< DataType, DIM> const * in_space,
             VectorSpace< DataType, DIM> const * out_space,
             const Vector& in_vec, 
             const Vector& out_vec)
  {
    assert (in_space != nullptr);
    assert (out_space != nullptr);
    
    std::vector<size_t> in_fe_inds;  
    std::vector<size_t> out_fe_inds;
    number_range<size_t>(0, 1, in_space->nb_fe(), in_fe_inds);
    number_range<size_t>(0, 1, out_space->nb_fe(), out_fe_inds);
    
    this->init(in_space, out_space, in_vec, out_vec, in_fe_inds, out_fe_inds);
  }
  
  void init (VectorSpace< DataType, DIM> const * in_space,
             VectorSpace< DataType, DIM> const * out_space,
             const Vector& in_vec, 
             const Vector& out_vec,
             const std::vector<size_t>& in_fe_inds, 
             const std::vector<size_t>& out_fe_inds)
  {
    assert (in_space != nullptr);
    assert (out_space != nullptr);
    
    NodalCellInterpolator * cell_inter = new NodalCellInterpolator(*out_space);
      
    this->init_with_linear_map_without_comm (cell_inter,
                                             in_space, 
                                             out_space,
                                             in_vec,
                                             out_vec, 
                                             in_fe_inds,
                                             out_fe_inds);
                                            
    delete cell_inter;
  }

  void interpolate (const Vector& in_vec, Vector& out_vec) const
  {
    //this->interpolate_with_linear_map_without_comm(in_vec, out_vec);
    this->interpolate_with_matrix_without_comm(in_vec, out_vec);
  }
  
  void interpolate_transpose (const Vector& in_vec, Vector& out_vec) const
  {
    this->interpolate_trans_with_matrix_without_comm(in_vec, out_vec);
  }
};

template < class LAD, int DIM> 
class FeInterMapRedNodal : public FeInterMapBase<LAD, DIM > 
{
  typedef typename LAD::DataType DataType;
  typedef typename LAD::VectorType Vector;
  typedef FeInterCellNodal<DataType, DIM, FeEvalLocal< DataType, DIM > >  NodalCellInterpolator;
  
public:

  FeInterMapRedNodal()
  : FeInterMapBase<LAD, DIM >()
  {
  }

  ~FeInterMapRedNodal()
  {
  }

  void init (VectorSpace< DataType, DIM> const * in_space,
             VectorSpace< DataType, DIM> const * out_space)
  {
    assert (in_space != nullptr);
    assert (out_space != nullptr);
    
    std::vector<size_t> in_fe_inds;  
    std::vector<size_t> out_fe_inds;
    number_range<size_t>(0, 1, in_space->nb_fe(), in_fe_inds);
    number_range<size_t>(0, 1, out_space->nb_fe(), out_fe_inds);
    
    this->init(in_space, out_space, in_fe_inds, out_fe_inds);
  }
  
  void init (VectorSpace< DataType, DIM> const * in_space,
             VectorSpace< DataType, DIM> const * out_space,
             const std::vector<size_t>& in_fe_inds, 
             const std::vector<size_t>& out_fe_inds)
  {
    assert (in_space != nullptr);
    assert (out_space != nullptr);
    
    this->init_without_linear_map_without_comm(in_space, 
                                               out_space, 
                                               in_fe_inds,
                                               out_fe_inds);
  }

  void interpolate (const Vector& in_vec, Vector& out_vec) const
  {
    NodalCellInterpolator * cell_inter = new NodalCellInterpolator(*this->out_space_);
    
    this->interpolate_without_linear_map_without_comm(cell_inter, in_vec, out_vec);
    delete cell_inter;
  }
  
};

////////////////////////////////////////////////////
///////////// FeInterMap ///////////////////////////
////////////////////////////////////////////////////

template < class LAD, int DIM >
FeInterMapBase<LAD, DIM>::FeInterMapBase()
: in_space_(nullptr), out_space_(nullptr), initialized_(false), parcom_(nullptr), 
  diag_(nullptr),
  diagT_(nullptr),
  odiag_(nullptr),
  odiagT_(nullptr)
{
}

template < class LAD, int DIM >
void FeInterMapBase<LAD, DIM>::clear()
{
  this->weights_.clear();
  this->in_dof_ids_.clear();
  this->out_dof_ids_.clear();
  this->cell_map_.clear();
  this->in_fe_inds_.clear();
  this->out_fe_inds_.clear();
  this->in_space_ = nullptr;
  this->out_space_ = nullptr;
  this->initialized_ = false;
  if (this->parcom_ != nullptr)
  {
    delete this->parcom_;
    this->parcom_ = nullptr;
  }
  if (this->diag_ != nullptr)
  {
    delete this->diag_;
    this->diag_ = nullptr;
  }
  if (this->odiag_ != nullptr)
  {
    delete this->odiag_;
    this->odiag_ = nullptr;
  }
  if (this->diagT_ != nullptr)
  {
    delete this->diagT_;
    this->diagT_ = nullptr;
  }
  if (this->odiagT_ != nullptr)
  {
    delete this->odiagT_;
    this->odiagT_ = nullptr;
  }
}

template < class LAD, int DIM >
template < class CellInterpolator >
void FeInterMapBase<LAD, DIM>::init_with_linear_map_without_comm (CellInterpolator * cell_inter,
                                                                  VectorSpace< DataType, DIM> const * in_space,
                                                                  VectorSpace< DataType, DIM> const * out_space,
                                                                  const Vector& in_vec, 
                                                                  const Vector& out_vec,
                                                                  const std::vector<size_t>& in_fe_inds, 
                                                                  const std::vector<size_t>& out_fe_inds)
{
  assert (cell_inter != nullptr);
  assert (in_space != nullptr);
  assert (out_space != nullptr);
  assert (in_fe_inds.size() > 0);
  assert (out_fe_inds.size() > 0);
  assert (in_fe_inds.size() == out_fe_inds.size());
  
#ifndef NDEBUG
  for (size_t l=0; l<in_fe_inds.size(); ++l)
  {
    const size_t in_fe_ind = in_fe_inds[l];
    const size_t out_fe_ind = out_fe_inds[l];
    assert(in_fe_ind >= 0);
    assert(out_fe_ind >= 0);
    assert(in_fe_ind < in_space->nb_fe());
    assert(out_fe_ind < out_space->nb_fe());
  }
#endif

  this->clear();
  this->parcom_ = new ParCom(in_space->get_mpi_comm());
   
  this->in_space_ = in_space;
  this->out_space_ = out_space;
  
  this->in_fe_inds_ = in_fe_inds;
  this->out_fe_inds_ = out_fe_inds;

  // create map to obtain for K of out_mesh all adjacent K' of in_mesh
  //LOG_INFO("init linear map", " create cell map ");
  mesh::find_adjacent_cells<DataType, DIM> (out_space->meshPtr(), in_space->meshPtr(), this->cell_map_);
  
  const int tdim = out_space->meshPtr()->tdim();
  const size_t out_num_cell = out_space->meshPtr()->num_entities(tdim);

  this->weights_.clear();
  this->weights_.resize(out_num_cell);
  this->out_dof_ids_.resize(out_num_cell);
  this->in_dof_ids_.resize(out_num_cell);

  std::vector< std::vector<DataType> > coeff;
      
  // loop over all cells in out_mesh
  //LOG_INFO("init linear map", " compute coefficients ");
  
  assert (in_space->meshPtr()->num_entities(DIM) == in_space->fe_manager().fe_tank_size());
  assert (out_space->meshPtr()->num_entities(DIM) == out_space->fe_manager().fe_tank_size());
    
  for (int out_index = 0; out_index < out_num_cell; ++out_index)
  {
    mesh::Entity out_cell (this->out_space_->meshPtr(), tdim, out_index);
   
    //std::cout << " out cell " << out_index << " : num cells " << this->cell_map_[out_index].size() << std::endl;

    // loop over all considered variables 
    for (size_t l=0; l<this->in_fe_inds_.size(); ++l)
    {
      const size_t in_fe_ind = this->in_fe_inds_[l];
      const size_t out_fe_ind = this->out_fe_inds_[l];
      const size_t nb_dofs_on_out_cell = this->out_space_->fe_manager().get_fe(out_index, out_fe_ind)->nb_dof_on_cell();
            
      // create basis evaluation object
      LOG_DEBUG(2, "   create basis eval object for " << this->cell_map_[out_index].size() << " cells ");
      BasisEvalLocal<DataType, DIM> basis_eval (*this->in_space_, in_fe_ind, this->cell_map_[out_index]);
     
      // pass evaluable object to cell interpolator
      cell_inter->set_function(&basis_eval);
        
      // evaluate cell interpolator
      LOG_DEBUG(2, "   evaluator cell interpolator ");
      cell_inter->compute_fe_coeff (&out_cell, out_fe_ind, coeff); 
        
      LOG_DEBUG(2, "   data handling ");
      this->weights_[out_index].push_back(coeff);

      // get global dof ids of in_space
      std::vector<DofId> in_gl_ids;
      basis_eval.get_basis_ids(in_gl_ids);
      this->in_dof_ids_[out_index].push_back(in_gl_ids);
      
      assert ( coeff[0].size() == in_gl_ids.size() );
      
      // get global dof ids of out cell
      std::vector< DofId > out_gl_ids (nb_dofs_on_out_cell);
      this->out_space_->get_dof_indices(out_fe_ind, out_index, &out_gl_ids);
      this->out_dof_ids_[out_index].push_back(out_gl_ids);
      
      assert ( coeff.size() == nb_dofs_on_out_cell);
      //std::cout << " out cell " << out_index << " : " << coeff.size() << " " << in_gl_ids.size() << std::endl;
      //log_2d_array(coeff, std::cout, 2);
    }
  }
  this->build_interpolation_matrix(in_vec, out_vec);
  this->initialized_ = true;
}

template < class LAD, int DIM >
void FeInterMapBase<LAD, DIM>::init_without_linear_map_without_comm (VectorSpace< DataType, DIM> const * in_space,
                                                                     VectorSpace< DataType, DIM> const * out_space,
                                                                     const std::vector<size_t>& in_fe_inds, 
                                                                     const std::vector<size_t>& out_fe_inds)
{
  assert (in_space != nullptr);
  assert (out_space != nullptr);
  assert (in_fe_inds.size() > 0);
  assert (out_fe_inds.size() > 0);
  assert (in_fe_inds.size() == out_fe_inds.size());
  
#ifndef NDEBUG
  for (size_t l=0; l<in_fe_inds.size(); ++l)
  {
    const size_t in_fe_ind = in_fe_inds[l];
    const size_t out_fe_ind = out_fe_inds[l];
    assert(in_fe_ind >= 0);
    assert(out_fe_ind >= 0);
    assert(in_fe_ind < in_space->nb_fe());
    assert(out_fe_ind < out_space->nb_fe());
  }
#endif

  this->clear();
  this->parcom_ = new ParCom(in_space->get_mpi_comm());
   
  this->in_space_ = in_space;
  this->out_space_ = out_space;
  
  this->in_fe_inds_ = in_fe_inds;
  this->out_fe_inds_ = out_fe_inds;

  // create map to obtain for K of out_mesh all adjacent K' of in_mesh
  mesh::find_adjacent_cells<DataType, DIM> (out_space->meshPtr(), in_space->meshPtr(), this->cell_map_);
     
  this->initialized_ = true;
}

template < class LAD, int DIM >
void FeInterMapBase<LAD, DIM>::build_interpolation_matrix(const Vector& in_vec, 
                                                          const Vector& out_vec)
{
  // interpolation weights below eps are considered to be zero
  const DataType eps = 1e-16;
  const int tdim = this->out_space_->meshPtr()->tdim();
  const size_t out_num_cell = this->out_space_->meshPtr()->num_entities(tdim);
  const size_t num_fe = this->in_fe_inds_.size();
 
  const int diag_nrow = out_vec.size_local();
  const int diag_ncol = in_vec.size_local();
  
  const int odiag_nrow = out_vec.size_local();
  const int odiag_ncol = in_vec.size_local_ghost();

  const int diagT_nrow = in_vec.size_local();
  const int diagT_ncol = out_vec.size_local();
  
  const int odiagT_nrow = in_vec.size_local();
  const int odiagT_ncol = out_vec.size_local_ghost();
    
  const int mat_diag_nnz_est = diag_nrow * 20;
  const int mat_odiag_nnz_est = odiag_nrow * 20;
  const int matT_diag_nnz_est = diagT_nrow * 20;
  const int matT_odiag_nnz_est = odiagT_nrow * 20;
  
  SortedArray<DofId> visited_i;
  
  std::vector<DofId> diag_coo_i;
  std::vector<DofId> diag_coo_j;
  std::vector<DataType> diag_coo_v;

  std::vector<DofId> diagT_coo_i;
  std::vector<DofId> diagT_coo_j;
  std::vector<DataType> diagT_coo_v;

  std::vector<DofId> odiag_coo_i;
  std::vector<DofId> odiag_coo_j;
  std::vector<DataType> odiag_coo_v;

  std::vector<DofId> odiagT_coo_i;
  std::vector<DofId> odiagT_coo_j;
  std::vector<DataType> odiagT_coo_v;

  // diag and offdiagonal part of interpolation matrix
  diag_coo_i.reserve(mat_diag_nnz_est);
  diag_coo_j.reserve(mat_diag_nnz_est);
  diag_coo_v.reserve(mat_diag_nnz_est);

  odiag_coo_i.reserve(mat_odiag_nnz_est);
  odiag_coo_j.reserve(mat_odiag_nnz_est);
  odiag_coo_v.reserve(mat_odiag_nnz_est);

  // diag and offdiagonal part of transposed interpolation matrix
  diagT_coo_i.reserve(matT_diag_nnz_est);
  diagT_coo_j.reserve(matT_diag_nnz_est);
  diagT_coo_v.reserve(matT_diag_nnz_est);

  odiagT_coo_i.reserve(matT_odiag_nnz_est);
  odiagT_coo_j.reserve(matT_odiag_nnz_est);
  odiagT_coo_v.reserve(matT_odiag_nnz_est);
  
  //SortedArray<DofID> out_ids;
  //SortedArray<DofId> in_ids;
  
  //out_ids.reserve(2*num_row);
  //in_ids.reserve(2*num_col);
  
  std::vector<DofId> tmp_loc_in_ids;
  std::vector<DofId> tmp_gh_in_ids;
  std::vector<DofId> tmp_loc_out_ids;
  std::vector<DofId> tmp_gh_out_ids;
      
  std::vector< DataType > dof_factors;
  
  //this->in_ids_.clear();
  //this->in_ids_.resize(num_col);
  //this->out_ids_.clear();
  //this->out_ids_.resize(num_row);
  
  // compute interpolation matrices in COO form
  
  // loop over all out cell
  for (size_t out_index = 0; out_index < out_num_cell; ++out_index)
  {
    dof_factors.clear();
    this->out_space_->dof().get_dof_factors_on_cell(out_index, dof_factors);

    // loop over all fe's
    for (size_t l=0; l<num_fe; ++l)
    {
      const size_t out_fe_ind   = this->out_fe_inds_[l];
      const size_t num_out_dofs = this->out_dof_ids_[out_index][l].size();
      const size_t num_in_dofs  = this->in_dof_ids_[out_index][l].size();

      // for indexing dof_factors
      size_t out_start_dof_on_cell = 0;
      for (size_t k=0; k<out_fe_ind; ++k)
      { 
        out_start_dof_on_cell += this->out_space_->dof().nb_dofs_on_cell(k, out_index);
      }
      
      // input: get local and ghost dof ids for current cell and fe 
      tmp_loc_in_ids.clear();
      tmp_loc_in_ids.reserve(num_in_dofs);
      tmp_gh_in_ids.clear();
      tmp_gh_in_ids.reserve(num_in_dofs);
      
      this->in_space_->global_2_local_and_ghost(this->in_dof_ids_[out_index][l],
                                                tmp_loc_in_ids,
                                                tmp_gh_in_ids);
                        
      // output: get local and ghost dof ids for current cell and fe
      tmp_loc_out_ids.clear();
      tmp_loc_out_ids.reserve(num_out_dofs);
      tmp_gh_out_ids.clear();
      tmp_gh_out_ids.reserve(num_out_dofs);
      
      this->out_space_->global_2_local_and_ghost(this->out_dof_ids_[out_index][l],
                                                 tmp_loc_out_ids,
                                                 tmp_gh_out_ids);
      
      // loop over output dofs
      for (size_t i=0; i!= num_out_dofs; ++i)
      {
        // get global out dof id
        const DofId gl_i = this->out_dof_ids_[out_index][l][i];
               
        // visit all out ids only once
        bool found_i = visited_i.find_insert(gl_i);
        if (found_i)
        {
          continue;
        }
        
        // get local and ghost out id
        const DofId loc_i = tmp_loc_out_ids[i];
        const DofId gh_i = tmp_gh_out_ids[i];
  
        assert (loc_i >= 0 || gh_i >= 0);
        const bool diag_i = (gh_i == -1);
          
        //this->out_ids_[loc_i] = gl_i;
        //out_ids.insert(loc_i);
        
        // compare VectorSpace::insert_dof_values()
        const DataType factor_i = 1. / dof_factors[out_start_dof_on_cell + i];
                
        // loop over in dofs
        for (size_t j=0; j!=num_in_dofs; ++j)
        {
          const DofId loc_j = tmp_loc_in_ids[j];
          const DofId gh_j = tmp_gh_in_ids[j];
          
          assert (loc_j >= 0 || gh_j >= 0);
          const bool diag_j = (gh_j == -1);
          
          DataType val = this->weights_[out_index][l][i][j] * factor_i;
          
          // skip entries close to zero
          if (std::abs(val) < eps)
          {
            continue;
          }
          
          // correct rounding errors if val is close to -1, 1
          if (std::abs(val-1.) < eps)
          {
            val = 1.;
          }
          else if (std::abs(val+1.) < eps)
          {
            val = -1.;
          }
          
          if (diag_i)
          {
            // out id is in diagonal part
            if (diag_j)
            {
              // in id is in diagonal part
              assert (loc_i < diag_nrow);
              assert (loc_j < diag_ncol);
              
              diag_coo_i.push_back(loc_i);
              diag_coo_j.push_back(loc_j);
              diag_coo_v.push_back(val);

              assert (loc_j < diagT_nrow);
              assert (loc_i < diagT_ncol);
              
              diagT_coo_i.push_back(loc_j);
              diagT_coo_j.push_back(loc_i);
              diagT_coo_v.push_back(val);
            }
            else
            {
              // in id is in offdiagonal part
              assert (loc_i < odiag_nrow);
              assert (gh_j < odiag_ncol);
              
              odiag_coo_i.push_back(loc_i);
              odiag_coo_j.push_back(gh_j);
              odiag_coo_v.push_back(val);
            }
          }
          else
          {
            // out id is in offdiagonal part
            if (diag_j)
            {
              // in id is in diagonal part
              assert (loc_j < odiagT_nrow);
              assert (gh_i < odiagT_ncol);
              
              odiagT_coo_i.push_back(loc_j);
              odiagT_coo_j.push_back(gh_i);
              odiagT_coo_v.push_back(val);
            }
          }
        }
      }
    }
  }
  
  const int diag_nnz   = diag_coo_i.size();
  const int odiag_nnz  = odiag_coo_i.size();
  const int diagT_nnz  = diagT_coo_i.size();
  const int odiagT_nnz = odiagT_coo_i.size();
  
  assert (diag_nnz > 0);
  assert (diagT_nnz > 0);
  
#ifdef nWITH_MKL
  this->diag_   = la::init_matrix<DataType>(diag_nnz,  diag_nrow,  diag_ncol,  "inter_diag",  la::CPU, la::MKL, la::CSR);
  this->diagT_  = la::init_matrix<DataType>(diagT_nnz, diagT_nrow, diagT_ncol, "inter_diagT", la::CPU, la::MKL, la::CSR);
  this->odiag_  = la::init_matrix<DataType>(odiag_nnz, odiag_nrow, odiag_ncol, "inter_odiag", la::CPU, la::MKL, la::CSR);
  this->odiagT_ = la::init_matrix<DataType>(odiagT_nnz,odiagT_nrow,odiagT_ncol,"inter_odiagT",la::CPU, la::MKL, la::CSR);
#else
  this->diag_   = la::init_matrix<DataType>(diag_nnz,  diag_nrow,  diag_ncol,  "inter_diag",  la::CPU, la::NAIVE, la::CSR);
  this->diagT_  = la::init_matrix<DataType>(diagT_nnz, diagT_nrow, diagT_ncol, "inter_diagT", la::CPU, la::NAIVE, la::CSR);
  this->odiag_  = la::init_matrix<DataType>(odiag_nnz, odiag_nrow, odiag_ncol, "inter_odiag", la::CPU, la::NAIVE, la::CSR);
  this->odiagT_ = la::init_matrix<DataType>(odiagT_nnz,odiagT_nrow,odiagT_ncol,"inter_odiagT",la::CPU, la::NAIVE, la::CSR);
#endif
  
  dynamic_cast<CPU_CSR_lMatrix <DataType> * >(this->diag_)  ->TransformFromCOO(&(diag_coo_i[0]),  &(diag_coo_j[0]),  &(diag_coo_v[0]),  diag_nrow,  diag_ncol,  diag_nnz); 
  dynamic_cast<CPU_CSR_lMatrix <DataType> * >(this->diagT_) ->TransformFromCOO(&(diagT_coo_i[0]), &(diagT_coo_j[0]), &(diagT_coo_v[0]), diagT_nrow, diagT_ncol, diagT_nnz); 
  dynamic_cast<CPU_CSR_lMatrix <DataType> * >(this->odiag_) ->TransformFromCOO(&(odiag_coo_i[0]), &(odiag_coo_j[0]), &(odiag_coo_v[0]), odiag_nrow, odiag_ncol, odiag_nnz); 
  dynamic_cast<CPU_CSR_lMatrix <DataType> * >(this->odiagT_)->TransformFromCOO(&(odiagT_coo_i[0]),&(odiagT_coo_j[0]),&(odiagT_coo_v[0]),odiagT_nrow,odiagT_ncol,odiagT_nnz); 
}

template < class LAD, int DIM >
void FeInterMapBase<LAD, DIM>::interpolate_with_matrix_without_comm (const Vector& in_vec, 
                                                                     Vector& out_vec) const
{
  assert (this->initialized_);  
  assert (this->diag_ != nullptr);
  assert (this->odiag_ != nullptr);
  assert (in_vec.interior().get_size() == this->diag_->get_num_col());
  assert (out_vec.interior().get_size() == this->diag_->get_num_row());
  assert (in_vec.ghost().get_size() == this->odiag_->get_num_col());
  assert (out_vec.interior().get_size() == this->odiag_->get_num_row());
  
  this->diag_->VectorMult(in_vec.interior(), &(out_vec.interior()));

  this->odiag_->VectorMultAdd(in_vec.ghost(), &(out_vec.interior()));
  
  out_vec.store_interior();
  out_vec.Update();
}

template < class LAD, int DIM >
void FeInterMapBase<LAD, DIM>::interpolate_trans_with_matrix_without_comm (const Vector& out_vec, 
                                                                           Vector& in_vec) const
{
  assert (this->initialized_);  
  assert (this->diagT_ != nullptr);
  assert (this->odiagT_ != nullptr);
  
  this->diagT_->VectorMult(out_vec.interior(), &(in_vec.interior()));

  this->odiagT_->VectorMultAdd(out_vec.ghost(), &(in_vec.interior()));

  in_vec.store_interior();
  in_vec.Update();
}

template < class LAD, int DIM >
void FeInterMapBase<LAD, DIM>::interpolate_with_linear_map_without_comm (const Vector& in_vec, 
                                                                         Vector& out_vec) const
{
  const int tdim = this->out_space_->meshPtr()->tdim();
  const size_t out_num_cell = this->out_space_->meshPtr()->num_entities(tdim);
  const size_t num_fe = this->in_fe_inds_.size();

  assert (this->initialized_);  
  assert (this->weights_.size() == out_num_cell);
  assert (this->out_dof_ids_.size() == out_num_cell);
  assert (this->in_dof_ids_.size() == out_num_cell);
  
  // loop over all cells in out_mesh
  for (size_t out_index = 0; out_index < out_num_cell; ++out_index)
  {
    assert (this->weights_[out_index].size() == num_fe);
    assert (this->in_dof_ids_[out_index].size() == num_fe);
    assert (this->out_dof_ids_[out_index].size() == num_fe);
          
    // loop over all considered variables 
    for (size_t l=0; l<num_fe; ++l)
    {
      const size_t out_fe_ind = this->out_fe_inds_[l];
      const size_t num_dofs_on_out_cell = this->out_dof_ids_[out_index][l].size();
      const size_t num_dofs_on_in_cell = this->in_dof_ids_[out_index][l].size();
        
      // get values from in_vector
      this->in_vals_.clear();
      this->in_vals_.resize(num_dofs_on_in_cell, 0.);
      in_vec.GetValues (&this->in_dof_ids_[out_index][l][0], num_dofs_on_in_cell, &in_vals_[0]);
      
      assert (in_vals_.size() == num_dofs_on_in_cell);
      assert (this->weights_[out_index][l].size() == num_dofs_on_out_cell);

      this->out_dof_values_.clear();
      this->out_dof_values_.resize(num_dofs_on_out_cell, 0.);

      // loop over all dofs on out_cell
      for (size_t i=0; i<num_dofs_on_out_cell; ++i)
      {
        const DofId gl_i = this->out_dof_ids_[out_index][l][i];
        
        assert (num_dofs_on_in_cell == this->weights_[out_index][l][i].size());
        
        // interpolate only locally owned dofs
        if (this->parcom_->rank() != this->out_space_->dof().owner_of_dof(gl_i))
        {
          continue;
        }
        
        DataType val_i = 0.;
                
        // loop over all dofs on in_cell
        for (size_t j=0; j<num_dofs_on_in_cell; ++j)
        {
          val_i += in_vals_[j] * this->weights_[out_index][l][i][j];
        }
        
        out_dof_values_[i] = val_i;
        
        //std::cout << out_index << " " << l << " : " << gl_i << " " << val_i << " " << in_vec.GetValue(gl_i) << std::endl;
        //out_vec.SetValue(gl_i, val_i);
      }
      this->out_space_->insert_dof_values(out_fe_ind, out_index, out_vec, out_dof_values_);
    }
  }
  out_vec.Update();
}

template < class LAD, int DIM >
template < class CellInterpolator >
void FeInterMapBase<LAD, DIM>::interpolate_without_linear_map_without_comm (CellInterpolator * cell_inter,
                                                                            const Vector& in_vec, 
                                                                            Vector& out_vec) const
{
  const int tdim = this->out_space_->meshPtr()->tdim();
  const size_t out_num_cell = this->out_space_->meshPtr()->num_entities(tdim);
  const size_t num_fe = this->in_fe_inds_.size();
  std::vector< std::vector<DataType> > cell_coeff;
  std::vector<DataType> dof_values;
  
  assert (this->initialized_);  

  // create Fe evaluator objects
  std::vector< FeEvalLocal<DataType, DIM>* > in_fe_evals (num_fe);
  for (size_t l=0; l<num_fe; ++l)
  {
    in_fe_evals[l] = new FeEvalLocal<DataType, DIM>(*this->in_space_, in_vec, this->in_fe_inds_[l]);
  }
  
  // loop over all cells in out_mesh
  for (int out_index = 0; out_index < out_num_cell; ++out_index)
  {         
    mesh::Entity out_cell = this->out_space_->meshPtr()->get_entity(tdim, out_index);
                            
    // loop over all considered variables 
    for (size_t l=0; l<num_fe; ++l)
    {
      const size_t in_fe_ind = this->in_fe_inds_[l];
      const size_t out_fe_ind = this->out_fe_inds_[l];

      // set trial cells for accelerating point search
      in_fe_evals[l]->set_trial_cells(this->cell_map_.at(out_index));

      // evaluate interpolation on current cell
      cell_inter->set_function(in_fe_evals[l]);
      cell_inter->compute_fe_coeff (&out_cell, out_fe_ind, cell_coeff);

      assert (cell_coeff.size() == this->out_space_->nb_dof_on_cell(out_fe_ind, out_index) );
      assert (cell_coeff[0].size() == 1);
    
      // insert values into vector sol
      dof_values.resize(cell_coeff.size(), 0.); 
      for (size_t i=0; i<dof_values.size(); ++i)
      {
        dof_values[i] = cell_coeff[i][0];
      }
      this->out_space_->insert_dof_values(out_fe_ind, out_index, out_vec, dof_values);
    }
  }
  out_vec.Update();
  
  for (size_t l=0; l<num_fe; ++l)
  {
    delete in_fe_evals[l];
  }
}

} // namespace hiflow
#endif
