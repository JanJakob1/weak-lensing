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

#include "fe_interpolation.h"

#include "dof/dof_partition_global.h"
#include "dof/dof_impl/dof_container_lagrange.h"
#include "fem/fe_manager.h"
#include "fem/fe_reference.h"
#include "fem/cell_trafo/cell_transformation.h"
#include "linear_algebra/block_matrix.h"
#include "mesh/geometric_search.h"
#include "mesh/mesh.h"
#include "mesh/iterator.h"
#include "mesh/entity.h"
#include "mesh/periodicity_tools.h"

/* (for InnerProductinterpolator)
#include "linear_solver/cg.h"
#include "linear_solver/gmres.h"
#include "linear_solver/preconditioner_ilupp.h"
* #ifdef WITH_HYPRE
#include "linear_solver/hypre_boomer_amg.h"
#endif
*/
#include "space/vector_space.h"
#include "space/element.h"


namespace hiflow {


template < class LAD, int DIM >
FEInterpolation< LAD, DIM >::FEInterpolation()
    : eps_(1e-12), nb_fe_(0), initialized_(false), from_space_(nullptr),
      to_space_(nullptr), tdim_(0), gdim_(0), rank_(-1), reuse_(true) {
  this->degrees_.resize(nb_fe_);
  this->is_cg_.resize(nb_fe_);
  this->initialized_ = false;
}

template < class LAD, int DIM >
void FEInterpolation< LAD, DIM >::init(VectorSpace< DataType, DIM > *from_space,
                                       VectorSpace< DataType, DIM > *to_space,
                                       const VectorType &from_vector,
                                       bool reuse) 
{
  assert(from_space != nullptr);
  assert(to_space != nullptr);

  this->clear();
  this->reuse_ = reuse;
  this->from_space_ = from_space;
  this->to_space_ = to_space;
  this->from_mesh_ = from_space->meshPtr();
  this->to_mesh_ = to_space->meshPtr();
  this->tdim_ = from_mesh_->tdim();
  this->gdim_ = to_mesh_->gdim();
  this->nb_fe_ = from_space->nb_fe();

  assert(this->from_mesh_ != nullptr);

  const MPI_Comm &comm = from_space->get_mpi_comm();

  MPI_Comm_rank(comm, &this->rank_);

  // check if patch interpolation is possible for given space
  bool valid_space = this->check_space();
  if (!valid_space) {
    LOG_DEBUG(0, " SpacePatchInterpolation: input space is not valid ! ");
    exit(-1);
  }

  // get fe degrees of input space
  this->degrees_.resize(this->nb_fe_, 0);
  this->is_cg_.resize(this->nb_fe_, true);

  for (int v = 0; v < this->nb_fe_; ++v) 
  {
    this->is_cg_[v] = from_space->fe_manager().is_cont(v);
    this->degrees_[v] = from_space->fe_manager().get_fe(0, v)->max_deg() * 2;
  }
}

template < class LAD, int DIM > 
void FEInterpolation< LAD, DIM >::clear() {
  this->nb_fe_ = 0;
  this->eps_ = 1e-12;
  this->degrees_.clear();
  this->is_cg_.clear();
  this->tdim_ = 0;
  this->gdim_ = 0;
  this->from_mesh_ = nullptr;
  this->from_space_ = nullptr;
  this->to_mesh_ = nullptr;
  this->to_space_ = nullptr;
  this->rank_ = 0;
  this->initialized_ = false;
  this->reuse_ = true;
}

///////////////////////////////////////////////////////
//////// FEInterpolationLagrange //////////////////////
///////////////////////////////////////////////////////

template < class LAD, int DIM >
FEInterpolationLagrange< LAD, DIM >::FEInterpolationLagrange()
    : FEInterpolation< LAD, DIM >() {}

template < class LAD, int DIM >
void FEInterpolationLagrange< LAD, DIM >::clear() 
{
  FEInterpolation< LAD, DIM >::clear();
  this->dof_weights_.clear();
}

template < class LAD, int DIM >
void FEInterpolationLagrange< LAD, DIM >::init( VectorSpace< DataType, DIM > *from_space, 
                                                VectorSpace< DataType, DIM > *to_space,
                                                const VectorType &from_vector, bool reuse) 
{
  FEInterpolation< LAD, DIM >::init(from_space, to_space, from_vector, reuse);
  this->create_dof_mapping();
  this->initialized_ = true;
}

template < class LAD, int DIM >
bool FEInterpolationLagrange< LAD, DIM >::check_space() const 
{
  for (size_t l=0; l<this->nb_fe_; ++l)
  {
    if(this->to_space_->fe_manager().get_fe(0,l)->type()!= doffem::FE_TYPE_LAGRANGE)
    {
      return false;
    }
    if(this->from_space_->fe_manager().get_fe(0,l)->type()!= doffem::FE_TYPE_LAGRANGE)
    {
      return false;
    }
    if (this->to_space_->fe_manager().get_fe(0,l)->nb_comp() != 1)
    {
      return false;
    }
    if (this->from_space_->fe_manager().get_fe(0,l)->nb_comp() != 1)
    {
      return false;
    }
  }
  return true;
}

template < class LAD, int DIM >
void FEInterpolationLagrange< LAD, DIM >::create_dof_mapping() 
{
  this->dof_weights_.resize(this->nb_fe_);

  // setup search object for from-mesh
  mesh::GeometricSearch<DataType, DIM> *search;

  if (this->from_mesh_->is_rectangular()) 
  {
    LOG_DEBUG(1, "Use RecGridGeometricSearch");
    search = new mesh::RecGridGeometricSearch<DataType, DIM>(this->from_mesh_);
  } 
  else 
  {
    LOG_DEBUG(1, "Use GridGeometricSearch");
    search = new mesh::GridGeometricSearch<DataType, DIM>(this->from_mesh_);
  }

  std::vector< mesh::MasterSlave > period = this->to_mesh_->get_period();
  LOG_DEBUG(1, "To_mesh: Period size = " << period.size());

  assert(this->to_mesh_->get_period().size() ==
         this->from_mesh_->get_period().size());

  // loop over all variables
  for (int v = 0; v < this->nb_fe_; ++v) 
  {
    LOG_DEBUG(1, "Var index " << v << " / " << this->nb_fe_ - 1);

    int num_dofs_to = this->to_space_->dof().nb_dofs_local(v);

    // loop over all cells in to-mesh
    int cell_counter = 0;
    for (mesh::EntityIterator it = this->to_mesh_->begin(this->tdim_), 
                              end_it = this->to_mesh_->end(this->tdim_);
                              it != end_it; ++it) 
    {
      LOG_DEBUG(3, "Cell index "
                       << cell_counter << " / "
                       << this->to_mesh_->num_entities(this->tdim_) - 1);

      // get dof indices for current cell w.r.t. to_space
      std::vector< int > global_dof_ids_on_cell_to;
      this->to_space_->get_dof_indices(v, *it, &global_dof_ids_on_cell_to);
      int num_dofs_on_cell_to = global_dof_ids_on_cell_to.size();
      std::vector< DataType > values;
      values.resize(num_dofs_on_cell_to, 0.);

      // get coordinates of Lagrange dof points on reference cell
      doffem::DofContainerLagrange<DataType, DIM> const * dofs_lagrange =
        dynamic_cast< doffem::DofContainerLagrange<DataType, DIM> const * > (this->to_space_->fe_manager().get_fe(it->index(),v)->dof_container());
      
      assert (dofs_lagrange != 0);
      std::vector< Coord > dof_coords_on_ref_cell = dofs_lagrange->get_dof_coords();
      assert (dof_coords_on_ref_cell.size() == num_dofs_on_cell_to);

      // map coords from reference cell to physical cell
      std::vector< Coord > dof_coords_on_phys_cell (dof_coords_on_ref_cell.size());
      std::vector< DataType > tmp_coords_to (dof_coords_on_phys_cell.size() * DIM);
      doffem::CellTransformation<DataType, DIM> * cell_trafo = this->to_space_->fe_manager().get_cell_transformation(it->index());

      for (size_t q=0; q<dof_coords_on_ref_cell.size(); ++q)
      {
        cell_trafo->transform(dof_coords_on_ref_cell[q], dof_coords_on_phys_cell[q]);
        for (size_t d =0 ; d<DIM; ++d)
        {
          tmp_coords_to[q*DIM+d] = dof_coords_on_phys_cell[q][d];
        }
      }
      
      // take into account periodic boundary
      std::vector< DataType > unperiod_coords_to;
      if (period.size() > 0)
      {
          unperiod_coords_to = mesh::unperiodify (tmp_coords_to, DIM, period);
      }
      else
      {
          unperiod_coords_to = tmp_coords_to;
      }

      std::vector< Coord > coords_to (tmp_coords_to.size() / DIM);
      for (int q=0; q<coords_to.size(); ++q)
      {
        for (int d=0; d<DIM; ++d)
        {
          coords_to[q][d] = unperiod_coords_to[q*DIM+d];
        }
      }

      SortedArray< int > trial_cells_from;

      // loop over to-dofs on current cell
      for (int j = 0; j < num_dofs_on_cell_to; ++j) 
      {
        int global_j = global_dof_ids_on_cell_to[j];

        // consider only locally owned dofs
        if (this->rank_ != this->to_space_->dof().owner_of_dof(global_j)) 
        {
          continue;
        }

        // check if dof has been already considered
        if (this->dof_weights_[v].find(global_j) != this->dof_weights_[v].end()) 
        {
          continue;
        }

        // insert empty map
        std::map< int, DataType > map_j;

        // search for cells in from_mesh which contain to-dof
        std::vector< int > cells_from;
//      std::vector< Coord > ref_points;

        search->find_cell(coords_to[j], trial_cells_from.data(), cells_from/*, ref_points*/);

        int num_cells_from = cells_from.size();
        LOG_DEBUG(2, "Dof " << global_j << ": "
                            << string_from_range(coords_to[j].begin(),
                                                 coords_to[j].end())
                            << ", found " << num_cells_from
                            << " cells in from mesh ");

#ifndef NDEBUG
        if (num_cells_from == 0) {
          LOG_DEBUG(1, "Dof " << global_j << ": "
                              << string_from_range(coords_to[j].begin(),
                                                   coords_to[j].end())
                              << ", found " << num_cells_from
                              << " cells in from mesh ");
        }
#endif
        // if no cell is found, continue with next to-dof
        if (cells_from.size() == 0) 
        {
          continue;
        }

        // augment trial cells by found cells
        for (int k = 0; k < cells_from.size(); ++k) 
        {
          trial_cells_from.find_insert(cells_from[k]);
        }
        // trial_cells_from.insert ( trial_cells_from.end ( ), cells_from.begin
        // ( ), cells_from.end ( ) ); if (print) std::cout << "dof_j " <<
        // global_j << " with coords " << coords_to[j][0] << ", " <<
        // coords_to[j][1] << std::endl;

        // loop over all shape functions of from-space with support in cells_j
        if (num_cells_from > 1) 
        {
          num_cells_from = 1;
        }

        for (int c = 0; c < num_cells_from; ++c) 
        {
          // get from-cell entity and fe-type on cell
          int cell_index_from = cells_from[c];
          mesh::Entity from_cell = this->from_mesh_->get_entity(this->tdim_, cell_index_from);
            
#if 0 //(old code)
          std::vector< doffem::RefElement< DataType, DIM > const * > fe_type_from;
          fe_type_from.push_back(this->from_space_->fe_manager().get_fe(cell_index_from, v));

          // setup evalshapefunctions object
          EvalShapeFunctions< DIM, DataType > shape_fun(fe_type_from);

          // get reference coordinates for to-dof in current from-cell
          Vec< DIM, DataType > ref_coord;

          // check if reference coordinates are already provided by the grid search
          /* this searchis still buggy ...
          if (ref_points.size() > c && ref_points[c].size() == DIM) {
            for (int d = 0; d < DIM; ++d) {
              ref_coord[d] = ref_points[c][d];
            }
          } else {
*         */
          if (true)
          {
            // compute reference coordinates by inverse cell transformation
            CellTransformation< DataType, DIM > *cell_trafo = this->from_space_->fe_manager().get_cell_transformation(cell_index_from);    
            cell_trafo->inverse(coords_to[j], ref_coord);

#ifndef NDEBUG
            DataType diff_x =
                std::abs(coords_to[j][0] - cell_trafo->x(ref_coord));
            LOG_DEBUG(4, "Residual cell trafo x " << diff_x);
            assert(diff_x < 1e-10);
            if (DIM >= 2) {
              DataType diff_y =
                  std::abs(coords_to[j][1] - cell_trafo->y(ref_coord));
              LOG_DEBUG(4, "Residual cell trafo y " << diff_y);
              assert(diff_y < 1e-10);
              if (DIM >= 3) {
                DataType diff_z =
                    std::abs(coords_to[j][2] - cell_trafo->z(ref_coord));
                LOG_DEBUG(4, "Residual cell trafo z " << diff_z);
                assert(diff_z < 1e-10);
              }
            }
#endif
          }
#endif
          // evaluate from space shape functions physical point coords_to[j]
          Element<DataType, DIM> from_elem (*this->from_space_, cell_index_from);

          std::vector< DataType> shape_vals (from_elem.nb_comp(v) * from_elem.nb_dof(v), 0.);
          from_elem.N_fe(coords_to[j], v, shape_vals);

          // get dof values on from space
          std::vector< int > global_dof_ids_on_cell_from;
          this->from_space_->get_dof_indices(v, from_cell, &global_dof_ids_on_cell_from);
          int num_dofs_on_cell_from = global_dof_ids_on_cell_from.size();

          assert(from_elem.nb_dof(v) == num_dofs_on_cell_from);

          for (int i = 0; i < num_dofs_on_cell_from; ++i) 
          {
            int global_i = global_dof_ids_on_cell_from[i];

            // store result in interpolator-mapping
            if (std::abs(shape_vals[i]) > this->eps_) 
            {
              // if (print) std::cout << "interpolating dof " << global_i << "
              // with shape value " << shape_vals[i] << std::endl;
              map_j.insert(std::pair< int, DataType >(global_i, shape_vals[from_elem.get_fe(v)->iv2ind(i,0)]));
            }
          }
        }
        this->dof_weights_[v].insert( std::pair< int, std::map< int, DataType > >(global_j, map_j));
      }
      cell_counter++;
    }
  }

  LOG_DEBUG(1, "DOF mapping finished");
  delete search;
}

template < class LAD, int DIM >
void FEInterpolationLagrange< LAD, DIM >::interpolate( const VectorType &from_vec, VectorType &to_vec) 
{
  if (!this->initialized_) 
  {
    if (this->to_space_ == nullptr || this->from_space_ == nullptr) 
    {
      std::cout << "FEInterpolation is not initialized ! " << std::endl;
      exit(-1);
    } 
    else 
    {
      this->init(this->from_space_, this->to_space_, from_vec, this->reuse_);
    }
  }

  to_vec.Zeros();
  /*
  int in_begin = from_vec.ownership_begin ( );
  int in_end = from_vec.ownership_end ( );
  int out_begin = to_vec.ownership_begin ( );
  int out_end = to_vec.ownership_end ( );

  std::cout << "  Ownership range in-vector:    " << in_begin << " - " << in_end
  << std::endl; std::cout << "  Ownership range out-vector:    " << out_begin <<
  " - " << out_end << std::endl;
   */

  // loop over all variables
  for (int v = 0; v < this->nb_fe_; ++v) 
  {
    // loop over all dofs for current variable
    for (typename std::map< int, std::map< int, DataType > >::const_iterator it_j = this->dof_weights_[v].begin();
         it_j != this->dof_weights_[v].end(); ++it_j) 
    {
      int j = it_j->first;
      if (this->rank_ != this->to_space_->dof().owner_of_dof(j)) 
      {
        continue;
      }

      // loop over interpolating dofs of from-space
      for (typename std::map< int, DataType >::const_iterator it_i = it_j->second.begin(); it_i != it_j->second.end(); ++it_i) 
      {
        std::vector< int > i(1);
        i[0] = it_i->first;
        DataType coeff = it_i->second;

        std::vector< DataType > in_val(1, 0.);

        LOG_DEBUG(2, "dof j = " << j << " interpolated by dof i = " << i[0]
                                << " with weight " << coeff);
#ifndef NDEBUG
        if (std::abs(coeff) < 1e-12) {
          LOG_DEBUG(1, "dof j = " << j << " interpolated by dof i = " << i[0]
                                  << " with weight " << coeff);
        }
#endif
        from_vec.GetValues(&i[0], 1, &in_val[0]);
        to_vec.Add(j, coeff * in_val[0]);
      }
    }
  }
  to_vec.Update();
  interpolate_constrained_vector< DataType, DIM >(*this->to_space_, to_vec);
  to_vec.Update();

  if (!this->reuse_) 
  {
    this->dof_weights_.clear();
    this->initialized_ = false;
  }
}

template class FEInterpolationLagrange< la::LADescriptorCoupledD, 2 >;
template class FEInterpolationLagrange< la::LADescriptorBlock< la::LADescriptorCoupledD >, 2 >;
template class FEInterpolationLagrange< la::LADescriptorCoupledD, 3 >;
template class FEInterpolationLagrange< la::LADescriptorBlock< la::LADescriptorCoupledD >, 3 >;

#ifdef WITH_HYPRE
template class FEInterpolationLagrange< la::LADescriptorHypreD, 2 >;
template class FEInterpolationLagrange< la::LADescriptorBlock< la::LADescriptorHypreD >, 2 >;
template class FEInterpolationLagrange< la::LADescriptorHypreD, 3 >;
template class FEInterpolationLagrange< la::LADescriptorBlock< la::LADescriptorHypreD >, 3 >;
#endif

#if 0 // TODO: DEBUGGING
    template<class LAD, int DIM>
    void FEInterpolationInnerProduct<LAD, DIM>::init ( VectorSpace<DataType>* from_space, VectorSpace<DataType>* to_space, const VectorType& from_vector, bool reuse )
    {
        std::cout << "FEInterpolationInnerProduct is not properly tested yet " << std::endl;
        exit(-1);

        FEInterpolation<LAD, DIM>::init ( from_space, to_space, from_vector, reuse );

        if ( this->comm_ != MPI_COMM_NULL )
        {
            MPI_Comm_free ( &this->comm_ );
        }

        assert ( to_space->get_mpi_comm ( ) != MPI_COMM_NULL );

        // determine nb. of processes
        int nb_procs;
        int info = MPI_Comm_size ( to_space->get_mpi_comm ( ), &nb_procs );
        assert ( info == MPI_SUCCESS );
        assert ( nb_procs > 0 );

        // retrieve my rank
        int my_rank;
        info = MPI_Comm_rank ( to_space->get_mpi_comm ( ), &my_rank );
        assert ( info == MPI_SUCCESS );
        assert ( my_rank >= 0 );
        assert ( my_rank < nb_procs );

        info = MPI_Comm_split ( to_space->get_mpi_comm ( ), 0, my_rank, &( this->comm_ ) );
        assert ( info == MPI_SUCCESS );

        this->from_vector_ = &from_vector;

        // setup assembler for from_space
        this->setup_asm();

        // setup Linear algebra structures
        this->setup_LA();

        // assemble right hand side
        this->assemble_matrix();

        // setup linear solver
        if (this->user_solver_ == nullptr)
        {
            this->setup_solver();
        }
        else
        {
            this->user_solver_->SetupOperator ( this->mass_matrix_ );
        }

        this->initialized_ = true;
    }

    template<class LAD, int DIM>
    bool FEInterpolationInnerProduct<LAD, DIM>::check_space ( ) const
    {
        // TODO

        return true;
    }

    template<class LAD, int DIM>
    void FEInterpolationInnerProduct<LAD, DIM>::setup_asm ( )
    {
        assert (this->from_space_ != nullptr);
        assert (this->from_vector_ != nullptr);
        assert (this->local_asm_ != nullptr);

        this->local_asm_->set_from_space ( this->from_space_, this->from_vector_ );
    }

    template<class LAD, int DIM>
    void FEInterpolationInnerProduct<LAD, DIM>::setup_LA ( )
    {
        assert (this->to_space_ != nullptr);
        int num_var = this->to_space_->get_nb_var();
        assert (this->to_coupling_vars_.size() == num_var);

        SparsityStructure sparsity;
        global_asm_.compute_sparsity_structure ( *this->to_space_, sparsity, &this->to_coupling_vars_ );

        // Initialize linear algebra structures
        this->couplings_.Init ( this->comm_ );
        std::vector<int> global_offsets;
        std::vector<int> ghost_dofs;
        std::vector<int> ghost_offsets;

        this->to_space_->get_la_couplings ( global_offsets, ghost_dofs, ghost_offsets );
        this->couplings_.InitializeCouplings ( global_offsets, ghost_dofs, ghost_offsets );

        // Initialize matrices and vectors
        this->mass_matrix_.Init ( this->comm_, this->couplings_ );
        this->mass_matrix_.InitStructure ( vec2ptr ( sparsity.diagonal_rows ),
                             vec2ptr ( sparsity.diagonal_cols ),
                             sparsity.diagonal_rows.size ( ),
                             vec2ptr ( sparsity.off_diagonal_rows ),
                             vec2ptr ( sparsity.off_diagonal_cols ),
                             sparsity.off_diagonal_rows.size ( ) );

        this->rhs_.Init ( this->comm_, this->couplings_ );
        this->rhs_.Zeros ( );
    }

    template<class LAD, int DIM>
    void FEInterpolationInnerProduct<LAD, DIM>::assemble_matrix ( )
    {
        this->mass_matrix_.Zeros();
        this->global_asm_.assemble_matrix ( *this->to_space_, boost::ref ( *this->local_asm_ ), this->mass_matrix_ );

        if ( !this->to_dirichlet_dofs_.empty ( ) )
        {
            this->mass_matrix_.diagonalize_rows ( vec2ptr ( this->to_dirichlet_dofs_ ), this->to_dirichlet_dofs_.size ( ), 1. );
        }
    }

    template<class LAD, int DIM>
    void FEInterpolationInnerProduct<LAD, DIM>::assemble_rhs ( )
    {
        assert (this->from_vector_ != nullptr);

        this->local_asm_->set_from_vector ( this->from_vector_ );
        this->rhs_.Zeros();
        this->global_asm_.assemble_vector ( *this->to_space_, boost::ref ( *this->local_asm_ ), this->rhs_ );

        if ( !this->to_dirichlet_dofs_.empty ( ) )
        {
            rhs_.SetValues ( vec2ptr ( this->to_dirichlet_dofs_ ), this->to_dirichlet_dofs_.size ( ), vec2ptr ( this->to_dirichlet_values_ ) );
        }
        rhs_.Update ( );
    }

    template<class LAD, int DIM>
    void FEInterpolationInnerProduct<LAD, DIM>::setup_solver ( )
    {
        assert (this->to_space_ != nullptr);
        int num_var = this->to_space_->get_nb_var();

        const int max_it = 1000;
        const int max_size = 100;
        const DataType abs_tol = 1e-20;
        const DataType rel_tol = 1e-12;

        this->cg_.InitControl(max_it, abs_tol, rel_tol, 1e6);
        this->cg_.SetPrintLevel(0);
        this->cg_.SetReuse(true);
        this->cg_.SetupOperator ( this->mass_matrix_ );

#if 0
        amg_.Clear();
        this->amg_.SetPreconditioningParameters();

        this->amg_.SetNumFunctions (num_var);
        this->amg_.SetCycleType (1);
        this->amg_.InitControl (1,0,0);
        this->amg_.SetCycleRelaxType (0, 1);
        this->amg_.SetCycleRelaxType (0, 2);
        this->amg_.SetRelaxWt (0.5);
        this->amg_.SetInterpType (4);
        this->amg_.SetStrongThreshold (0.55);
        this->amg_.SetAggNumLevels (0);
        this->amg_.SetCoarsenType (6);
        this->amg_.SetCycleNumSweeps (1, 1);
        this->amg_.SetCycleNumSweeps (1, 2);
        this->amg_.SetSmoothType (0);
        this->amg_.SetSmoothNumLevels (-1);
        this->amg_.SetMaxCoarseSize (5);
        this->amg_.SetMaxLevels (10);
        this->amg_.SetCycleRelaxType (9, 3);
        this->amg_.SetCycleNumSweeps (1, 3);
        this->amg_.SetVariant (3);
        this->amg_.SetOverlap (1);
        this->amg_.SetDomainType (1);
        this->amg_.SetSchwarzUseNonSymm (0);

        std::vector<int> vars (num_var, 0);
        for (int d=0; d<num_var; ++d)
        {
            vars[d] = d;
        }
        this->to_dof_func_ = this->to_space_->get_dof_func ( vars );
        this->amg_.SetDofFunc ( this->dof_func_ );

        this->amg_.SetReuse(true);
        this->amg_.Init();
        this->amg_.SetupOperator ( this->mass_matrix_ );

        this->gmres_.InitParameter(max_size, "RightPreconditioning");
        this->gmres_.SetupPreconditioner (this->amg_ );
#else
#if 0
        const int prepro_type   = 0;
        const int precond_no    = 11;
        const int max_levels    = 20;
        const DataType mem_factor = 0.8;
        const DataType threshold  = 1.75;
        const DataType min_pivot  = 0.01;
        this->ilupp_.InitParameter(prepro_type, precond_no, max_levels, mem_factor, threshold, min_pivot);
        this->ilupp_.SetupOperator ( this->mass_matrix_ );

        this->gmres_.InitParameter(max_size, "RightPreconditioning");
        this->gmres_.SetupPreconditioner (this->ilupp_ );

#else
        this->cg_.InitParameter("NoPreconditioning");
#endif
#endif
    }


    template<class LAD, int DIM>
    void FEInterpolationInnerProduct<LAD, DIM>::interpolate ( const VectorType& from_vec, VectorType& to_vec )
    {
        std::cout << "FEInterpolationInnerProduct is not properly tested yet " << std::endl;
        exit(-1);

        if ( !this->initialized_ )
        {
            if (this->to_space_ == nullptr || this->from_space_ == nullptr)
            {
                std::cout << "FEInterpolation is not initialized ! " << std::endl;
                exit ( -1 );
            }
            else
            {
                this->init(this->from_space_, this->to_space_, from_vec, this->reuse_);
            }
        }

        // assemble rhs
        this->from_vector_ = &from_vec;
        this->assemble_rhs();

        // solve linear system
        to_vec.Zeros();
        this->cg_.Solve ( this->rhs_, &to_vec );
        to_vec.Update ( );
        interpolate_constrained_vector<LAD> ( *this->to_space_, to_vec );
        to_vec.Update ( );

        LOG_DEBUG(1, "CG ended after " << cg_.iter() << " iterations with residual norm " << cg_.res() );

        if (!this->reuse_)
        {
            this->mass_matrix_.Clear();
            this->rhs_.Clear();
            this->cg_.Clear();
            //this->amg_.Clear();
            //this->ilupp_.Clear();
            this->to_dof_func_.clear();
            this->initialized_ = false;
        }
    }

    template class FEInterpolationInnerProduct<LADescriptorCoupledD, 1>;
    template class FEInterpolationInnerProduct<LADescriptorCoupledD, 2>;
    template class FEInterpolationInnerProduct<LADescriptorCoupledD, 3>;

    template class FEInterpolationInnerProduct<LADescriptorBlock<LADescriptorCoupledD>, 1>;
    template class FEInterpolationInnerProduct<LADescriptorBlock<LADescriptorCoupledD>, 2>;
    template class FEInterpolationInnerProduct<LADescriptorBlock<LADescriptorCoupledD>, 3>;
#ifdef WITH_HYPRE
    template class FEInterpolationInnerProduct<LADescriptorHypreD, 1>;
    template class FEInterpolationInnerProduct<LADescriptorHypreD, 2>;
    template class FEInterpolationInnerProduct<LADescriptorHypreD, 3>;
    
    template class FEInterpolationInnerProduct<LADescriptorBlock<LADescriptorHypreD>, 1>;
    template class FEInterpolationInnerProduct<LADescriptorBlock<LADescriptorHypreD>, 2>;
    template class FEInterpolationInnerProduct<LADescriptorBlock<LADescriptorHypreD>, 3>;
#endif
#endif
} // namespace hiflow
