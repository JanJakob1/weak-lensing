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

/// \author Aksel Alpay, Martin Wlotzka

#ifndef GMG_LEVEL_H
#define GMG_LEVEL_H

#include <time.h>

#include "boost/lexical_cast.hpp"
#include "boost/shared_ptr.hpp"
#include <mpi.h>
#include <string>
#include <vector>
#include <cassert>
#include <stdexcept>

#include "assembly/global_assembler.h"
#include "assembly/standard_assembly.h"
#include "linear_algebra/la_descriptor.h"
#include "linear_algebra/lmp/lmatrix_csr_cpu.h" // for CPU_CSR_lMatrix
#include "linear_algebra/lmp/lvector_cpu.h"     // for CPU_lVector
#include "linear_solver/gpu-based/async_iter_gpu.h"
#include "linear_solver/linear_solver.h"
#include "mesh/mesh_tools.h"
#include "mesh/communication.h" // for SharedVertexTable
#include "visualization/cell_visualization.h"
#include "visualization/vtk_writer.h"
#include "space/fe_evaluation.h"
#include "tools/mpi_tools.h" // for mpi_data_type<DataType>::get_type()

#include "util.h"

#ifdef WITH_CUDA
#include "linear_algebra/lmp/cuda/cuda_gmg.h"
#endif

#include <csignal>

#ifndef NDEBUG
#include <fstream>
#endif

using namespace hiflow;
using namespace hiflow::la;
using namespace hiflow::mesh;
using namespace hiflow::doffem;

namespace hiflow {
namespace la {
namespace gmg {
/// Contains basic linear algebra settings

struct Settings {
  PLATFORM platform;
  IMPLEMENTATION implementation;
  MATRIX_FORMAT matrix_format;
};

/// This class represents a single level within the MultiLevelHierarchy.
/// It initializes a mesh, a vectorspace, matrices and vectors
/// and provides functions to assemble and solve the system on this level.
/// Note that a BasicSingleLevel object is rather a handle than an actual,
/// copyable object. Although the object can be copied, it will still
/// use the same pointers to the underlying, non-copying HiFlow data structures.
/// A deep copy of a BasicSingleLevel object is not possible.

template < class LAD, int DIM > class BasicLevel {
public:
  TYPE_FROM_CLASS(LAD, MatrixType);
  TYPE_FROM_CLASS(LAD, VectorType);
  TYPE_FROM_CLASS(LAD, DataType);

  typedef boost::shared_ptr< MatrixType > MatrixPtr;
  typedef boost::shared_ptr< VectorType > VectorPtr;
  typedef boost::shared_ptr< VectorSpace< DataType, DIM > > SpacePtr;
  typedef boost::shared_ptr< LaCouplings > CouplingsPtr;

  typedef StandardGlobalAssembler< DataType, DIM > GlobalAssemblerType;

  /// Initializes this level. Collective on the supplied communicator.
  /// @param global_comm The communicator that shall be used by this level
  /// @param master_mesh The master mesh that is to be distributed. It
  /// only needs to be a valid mesh pointer on the master rank. The
  /// BasicSingleLevel will then take care of its partitioning
  /// and distribution.
  /// This will be the coarsest mesh of the hierarchy. The MultiLevelHierarchy
  /// will obtain meshes for other levels of the hierarchy by refining this
  /// mesh.
  /// @param fe_degrees Specifies the fe-degrees. The i-th entry of the
  /// vector is the polynomial degree of the i-th variable of the problem
  /// @param global_asm A pointer to a global assembler object
  /// @param Settings An instance of the settings struct that provides
  /// information about the linear algebra routines and data structures
  /// to be used.

  BasicLevel(MPI_Comm global_comm, MPI_Comm partial_comm,
             const MeshPtr &local_mesh, int master_rank,
             const std::vector< int > &fe_degrees,
             StandardGlobalAssembler< DataType, DIM > *global_asm,
             const Settings &settings);

  virtual ~BasicLevel() {}

  /// @return The process rank of the calling process in the communicator
  /// of this level

  int rank(void) const { return global_rank_; }

  int partial_rank(void) const { return partial_rank_; }

  /// @return The communicator that is used by this level, or MPI_COMM_NULL
  /// if the calling process is contained in the communicator for this level

  MPI_Comm comm(void) const { return global_comm_; }

  MPI_Comm partial_comm(void) const { return partial_comm_; }

  /// @return A pointer to the mesh used by this level, or a nullptr pointer
  /// if the calling process is not contained in the communicator for this
  /// level.

  MeshPtr mesh(void) const { return mesh_; }

  /// @return A pointer to the vector space used by this level, or a nullptr
  /// pointer if the calling process is not contained in the communicator for
  /// this level.

  SpacePtr space(void) const { return space_; }

  CouplingsPtr couplings(void) const { return couplings_; }

  /// @return The solution vector of this level, or a nullptr pointer
  /// if the calling process is not contained in the communicator for this
  /// level.

  VectorPtr sol(void) const { return sol_; }

  /// @return The right hand side vector of this level, or a nullptr pointer
  /// if the calling process is not contained in the communicator for this
  /// level.

  VectorPtr rhs(void) const { return rhs_; }

  /// @return The matrix of this level, or a nullptr pointer if the calling process
  /// is not contained in the communicator for this level.

  MatrixPtr matrix(void) const { return matrix_; }

  MatrixPtr restriction_matrix(void) const { return restriction_matrix_; }

  MatrixPtr interpolation_matrix(void) const { return interpolation_matrix_; }

  VectorPtr tmp_vector(void) const { return tmp_; }

  /// Assembles the Matrix and right hand side vector of this level.
  /// Collective on this level's communicator.
  /// @param v_asm The vector assembling function
  /// @param m_asm The matrix assembling function
  void
  assemble_system(typename GlobalAssemblerType::VectorAssemblyFunction v_asm,
                  typename GlobalAssemblerType::MatrixAssemblyFunction m_asm);

  /// Assembles only the matrix of this level. Collective on this level's
  /// communicator
  /// @param m_asm The matrix assembling function
  void
  assemble_matrix(typename GlobalAssemblerType::MatrixAssemblyFunction m_asm);

  /// @return whether the calling process is contained in this level's
  /// communicator, i.e. whether it has been scheduled to do work on this
  /// level.

  inline bool is_scheduled_to_this_process() const {
    return global_comm_ != MPI_COMM_NULL;
  }

  /// Solves the system on this level. Collective on this level's communicator.
  /// @param solver The linear solver that shall be used.

  void solve_system(LinearSolver< LAD > &solver) {
    if (is_scheduled_to_this_process()) {
      solver.SetupOperator(*matrix_);
      solver.Solve(*rhs_, sol_.get());
    }
  }

  /// Writes a visualization of a vector of this level to a file. Collective
  /// on this level's communicator.
  /// The \c num_intervals parameter of the visualization (which determines
  /// how many cells in the visualization correspond to real cells) will
  /// be set to the maximum fe degree of all variables.
  /// @param v A pointer to a vector living on this level.
  /// @param filename The name of the file where the visualization is
  /// to be saved.
  void visualize_vector(const VectorPtr v, const std::string &filename) const;

  /// Writes a visualization of the solution vector to a file. Collective
  /// on this level's communicator.
  /// The \c num_intervals parameter of the visualization (which determines
  /// how many cells in the visualization correspond to real cells) will
  /// be set to the maximum fe degree of all variables.
  /// @param filename The name of the file where the visualization is
  /// to be saved.

  void visualize_solution(const std::string &filename) const {
    visualize_vector(sol(), filename);
  }

  /// Sets the solution to zero again.
  /// Collective on this level's communicator.

  void reset_solution(void) {
    if (is_scheduled_to_this_process())
      sol_->Zeros();
  }

  /// Sets the matrix, rhs vector and the solution to zero again.
  /// After a call to this function, assemble_system() must be executed
  /// to make the system solvable.
  /// Collective on this level's communicator.

  void reset_system(void) {
    if (is_scheduled_to_this_process()) {
      rhs_->Zeros();
      matrix_->Zeros();
      reset_solution();
    }
  }

  /// Creates a new vector using the vectorspace and settings of this level.

  VectorPtr create_vector(void) const;

  MatrixPtr create_matrix(void) const;

  MatrixPtr create_matrix(const LaCouplings &cp,
                          const SparsityStructure &sp) const;

  void create_restriction_matrix(void);

  void create_interpolation_matrix(void);

  /// Initializes Dirichlet boundary conditions
  /// @param eval The Dirichlet Evaluator
  /// @param correct_only_matrix whether dirichlet bc should only be applied to
  /// the matrix or to the right hand side as well.

  template < class DirichletEvaluator >
  void init_boundary_conditions(DirichletEvaluator &eval,
                                bool correct_only_matrix = false) {
    if (is_scheduled_to_this_process()) {
      dirichlet_dofs_.clear();
      dirichlet_values_.clear();

      for (int var = 0; var < space_->aaa_get_nb_var(); ++var)
        compute_dirichlet_dofs_and_values(eval, *space_, var, dirichlet_dofs_,
                                          dirichlet_values_);

      if (!dirichlet_dofs_.empty()) {
        if (!correct_only_matrix) {
          rhs_->SetValues(vec2ptr(dirichlet_dofs_), dirichlet_dofs_.size(),
                          vec2ptr(dirichlet_values_));
          sol_->SetValues(vec2ptr(dirichlet_dofs_), dirichlet_dofs_.size(),
                          vec2ptr(dirichlet_values_));
        }

        // Correct Dirichlet dofs.
        matrix_->diagonalize_rows(vec2ptr(dirichlet_dofs_),
                                  dirichlet_dofs_.size(), 1.0);
      }
    }
  }

  const std::vector< int > &dirichlet_dofs(void) const {
    return dirichlet_dofs_;
  }

  const std::vector< int > &pids(void) const { return pids_; }

protected:
  std::string class_name;

  /// Initializes the linear algebra part of this level
  void init_la(void);

  /// Initilizes the mesh of this level. Collective on this level's
  /// communicator. The mesh will be partitioned using METIS if WITH_METIS has
  /// been defined. Otherwise, the naive partitioner of HiFlow will be used.
  /// @param master_mesh The master mesh that shall be distributed. Only
  /// needs to be a valid argument on the master process.
  void init_mesh(const MeshPtr &master_mesh);

  /// function to be registered as signal handler for SIGUSR1, does nothing;
  /// default signal handler might terminate the process, but it shall just wake up

  MPI_Comm global_comm_;
  MPI_Comm partial_comm_;
  MeshPtr mesh_;
  SpacePtr space_;
  CouplingsPtr couplings_;
  VectorPtr rhs_, sol_, tmp_;
  MatrixPtr matrix_, restriction_matrix_, interpolation_matrix_;
  StandardGlobalAssembler< DataType, DIM > *global_asm_;

  std::vector< int > dirichlet_dofs_;
  std::vector< DataType > dirichlet_values_;

  std::vector< int > fe_degrees_;

  int master_rank_;

  Settings settings_;

  int global_size_;
  int global_rank_;
  int partial_rank_;

  std::vector< int > pids_;
};

/// A BasicLevel implementation that contains two pointers
/// to InterGridConnection objects to the next coarser and next finer grids.

template < class LAD, int DIM, class ConnectionType >
class ConnectedLevel : public BasicLevel< LAD, DIM > {
public:
  TYPE_FROM_CLASS(LAD, DataType);

  typedef ConnectionType Connection;
  /// Initializes this level without setting the pointers to the
  /// InterGridConnection objects.
  /// Collective on the supplied communicator.
  /// @param comm The communicator that shall be used by this level
  /// @param master_mesh The master mesh that is to be distributed. It
  /// only needs to be a valid mesh pointer on the master rank. The
  /// \c BasicSingleLevel will then take care of its partitioning
  /// and distribution.
  /// This will be the coarsest mesh of the hierarchy. The MultiLevelHierarchy
  /// will obtain meshes for other levels of the hierarchy by refining this
  /// mesh.
  /// @param fe_degrees Specifies the fe-degrees. The i-th entry of the
  /// vector is the polynomial degree of the i-th variable of the problem
  /// @param global_asm A pointer to a global assembler object
  /// @param Settings An instance of the settings struct that provides
  /// information about the linear algebra routines and data structures
  /// to be used.

  ConnectedLevel(MPI_Comm global_comm, MPI_Comm partial_comm,
                 const MeshPtr &master_mesh, unsigned master_rank,
                 const std::vector< int > &fe_degrees,
                 StandardGlobalAssembler< DataType, DIM > *global_asm,
                 const Settings &settings)
      : BasicLevel< LAD, DIM >(global_comm, partial_comm, master_mesh, master_rank,
                          fe_degrees, global_asm, settings) {}

  /// Initializes this level from a BasicSingleLevel object and pointers
  /// to InterGridConnection objects to the finer and coarser levels.
  /// @param lvl A BasicSingleLevel object from which the ConnectedSingleLevel
  /// will be initialized.
  /// Note that the object created by this constructor will not be
  /// independent of lvl, as the vectors, matrices, vectorspace and mesh
  /// pointers will be shared.
  /// These are non-copyable objects, therefore a deep copy is not possible
  /// (or at least not worth the effort)
  /// @param to_finer A pointer to an InterGridConnection object linking
  /// this level with the next finer level.
  /// @param to_coarser A pointer to an InterGridConnection object linking
  /// this level with the next coarser level.

  ConnectedLevel(const BasicLevel< LAD, DIM > &lvl,
                 const boost::shared_ptr< Connection > &to_finer,
                 const boost::shared_ptr< Connection > &to_coarser)
      : BasicLevel< LAD, DIM >(lvl), connection_to_finer_grid_(to_finer),
        connection_to_coarser_grid_(to_coarser) {}

  virtual ~ConnectedLevel() {}

  /// @return A pointer to the InterGridConnection object connecting
  /// this level with the next finer grid.

  boost::shared_ptr< ConnectionType >
  get_connection_to_next_finer_grid() const {
    return connection_to_finer_grid_;
  }

  /// @return A pointer to the InterGridConnection object connecting
  /// this level with the next coarser grid.

  boost::shared_ptr< ConnectionType >
  get_connection_to_next_coarser_grid() const {
    return connection_to_coarser_grid_;
  }

  /// Sets the connections of this level to the next coarser and finer levels.
  /// @param to_coarser A pointer to an InterGridConnection object connecting
  /// this level with the next coarser level.
  /// @param to_finer A pointer to an InterGridConnection object connecting
  /// this level with the next finer level.

  void set_connections(const boost::shared_ptr< ConnectionType > &to_coarser,
                       const boost::shared_ptr< ConnectionType > &to_finer) {
    connection_to_coarser_grid_ = to_coarser;
    connection_to_finer_grid_ = to_finer;
  }

private:
  boost::shared_ptr< ConnectionType > connection_to_finer_grid_;
  boost::shared_ptr< ConnectionType > connection_to_coarser_grid_;
};

/// A connected \c ConnectedLevel that additionally contains a vector to store
/// the residual, as required by the geometric multigrid algorithm.

template < class LAD, int DIM, class ConnectionType >
class GMGLevel : public ConnectedLevel< LAD, DIM, ConnectionType > {
public:
  TYPE_FROM_CLASS(LAD, DataType);

  typedef ConnectedLevel< LAD, DIM, ConnectionType > BaseType;
  IMPORT_FROM_BASECLASS(BaseType, VectorPtr);

  //   using ConnectedLevel<LAD, ConnectionType>::create_vector;

  /// Initializes this level without setting the pointers to the
  /// InterGridConnection objects.
  /// Collective on the supplied communicator.
  /// @param global_comm The communicator that shall be used by this level
  /// @param master_mesh The master mesh that is to be distributed. It
  /// only needs to be a valid mesh pointer on the master rank. The
  /// \c BasicSingleLevel will then take care of its partitioning
  /// and distribution.
  /// This will be the coarsest mesh of the hierarchy. The MultiLevelHierarchy
  /// will obtain meshes for other levels of the hierarchy by refining this
  /// mesh.
  /// @param fe_degrees Specifies the fe-degrees. The i-th entry of the
  /// vector is the polynomial degree of the i-th variable of the problem
  /// @param global_asm A pointer to a global assembler object
  /// @param Settings An instance of the settings struct that provides
  /// information about the linear algebra routines and data structures
  /// to be used.

  GMGLevel(MPI_Comm global_comm, MPI_Comm partial_comm,
           const MeshPtr &master_mesh, unsigned master_rank,
           const std::vector< int > &fe_degrees,
           StandardGlobalAssembler< DataType, DIM > *global_asm,
           const Settings &settings)
      : ConnectedLevel< LAD, DIM, ConnectionType >(global_comm, partial_comm,
                                              master_mesh, master_rank,
                                              fe_degrees, global_asm, settings)
#ifdef WITH_CUDA
        ,
        R_diag_(0), R_offdiag_(0), A_diag_(0), A_offdiag_(0), P_diag_(0),
        P_offdiag_(0), lvec_sol_(0), lvec_sol_ghost_(0), lvec_rhs_(0),
        lvec_rhs_ghost_(0), lvec_res_(0), lvec_res_ghost_(0), dev_sol_(0),
        dev_sol_ghost_(0), dev_rhs_(0), dev_rhs_ghost_(0), dev_res_(0),
        dev_res_ghost_(0), dev_tmp_(0), host_sol_(0), host_sol_ghost_(0),
        host_rhs_(0), host_rhs_ghost_(0), host_res_(0), host_res_ghost_(0),
        row_R_diag_(0), row_R_offdiag_(0), row_A_diag_(0), row_A_offdiag_(0),
        row_P_diag_(0), row_P_offdiag_(0), col_R_diag_(0), col_R_offdiag_(0),
        col_A_diag_(0), col_A_offdiag_(0), col_P_diag_(0), col_P_offdiag_(0),
        nnz_A_diag_(0), nnz_A_offdiag_(0), nnz_R_diag_(0), nnz_R_offdiag_(0),
        nnz_P_diag_(0), nnz_P_offdiag_(0), nrows_(0), block_dim_(128),
        grid_transfer_on_device_(false), async_iter_gpu_(0)
#endif
  {
    res_ = this->create_vector();
  }

  virtual ~GMGLevel() { clear(); }

  /// @return The residual vector of this level

  VectorPtr res(void) const { return res_; }

  /// Computes the restricted residual from the current solution,
  /// i.e. r = R(b-Ax).
  void prepare_grid_transfer_on_device(LinearSolver< LAD > *smoother);

  void CopySolToDevice() {
    if (this->is_scheduled_to_this_process()) {
      assert(grid_transfer_on_device_);
#ifdef WITH_CUDA
      if (async_iter_gpu_ == 0) {
        memcpy(host_sol_, lvec_sol_->buffer, nrows_ * sizeof(DataType));
        cudaSetDevice(0);
        cudaMemcpyAsync(dev_sol_, host_sol_, nrows_ * sizeof(DataType),
                        cudaMemcpyHostToDevice, stream_);
      }
#endif
    }
  }

  void CopySolGhostToDevice(void) {
    if (this->is_scheduled_to_this_process()) {
      assert(grid_transfer_on_device_);
#ifdef WITH_CUDA
      if (nnz_A_offdiag_ > 0) {
        memcpy(host_sol_ghost_, lvec_sol_ghost_->buffer,
               lvec_sol_ghost_->get_size() * sizeof(DataType));
        cudaSetDevice(0);
        cudaMemcpyAsync(dev_sol_ghost_, host_sol_ghost_,
                        lvec_sol_ghost_->get_size() * sizeof(DataType),
                        cudaMemcpyHostToDevice, stream_);
      }
#endif
    }
  }

  void CopyRhsToDevice() {
    if (this->is_scheduled_to_this_process()) {
      assert(grid_transfer_on_device_);
#ifdef WITH_CUDA
      if (async_iter_gpu_ == 0) {
        memcpy(host_rhs_, lvec_rhs_->buffer, nrows_ * sizeof(DataType));
        cudaSetDevice(0);
        cudaMemcpyAsync(dev_rhs_, host_rhs_, nrows_ * sizeof(DataType),
                        cudaMemcpyHostToDevice, stream_);
      }
#endif
    }
  }

  void CopyResToDevice() {
    if (this->is_scheduled_to_this_process()) {
      assert(grid_transfer_on_device_);
#ifdef WITH_CUDA
      memcpy(host_res_, lvec_res_->buffer, nrows_ * sizeof(DataType));
      cudaSetDevice(0);
      cudaMemcpyAsync(dev_res_, host_res_, nrows_ * sizeof(DataType),
                      cudaMemcpyHostToDevice, stream_);
#endif
    }
  }

  void CopyResGhostToDevice(void) {
    if (this->is_scheduled_to_this_process()) {
      assert(grid_transfer_on_device_);
#ifdef WITH_CUDA
      if ((nnz_P_offdiag_ > 0) || (nnz_R_offdiag_ > 0)) {
        memcpy(host_res_ghost_, lvec_res_ghost_->buffer,
               lvec_res_ghost_->get_size() * sizeof(DataType));
        cudaSetDevice(0);
        cudaMemcpyAsync(dev_res_ghost_, host_res_ghost_,
                        lvec_res_ghost_->get_size() * sizeof(DataType),
                        cudaMemcpyHostToDevice, stream_);
      }
#endif
    }
  }

  void ComputeResidualOnDevice(void) {
    if (this->is_scheduled_to_this_process()) {
      assert(grid_transfer_on_device_);
#ifdef WITH_CUDA
      cudaSetDevice(0);

      if (nnz_A_offdiag_ > 0) {
        cuda_gmg_compute_residual(A_diag_, col_A_diag_, row_A_diag_, nrows_,
                                  A_offdiag_, col_A_offdiag_, row_A_offdiag_,
                                  dev_sol_, dev_rhs_, dev_res_, dev_sol_ghost_,
                                  grid_dim_, block_dim_, stream_);
      } else {
        cuda_gmg_compute_residual_no_offdiag(
            A_diag_, col_A_diag_, row_A_diag_, nrows_, dev_sol_, dev_rhs_,
            dev_res_, grid_dim_, block_dim_, stream_);
      }
#endif
    }
  }

  void ComputeRestrictionOnDevice(void) {
    if (this->is_scheduled_to_this_process()) {
      assert(grid_transfer_on_device_);
#ifdef WITH_CUDA
      cudaSetDevice(0);

      if (nnz_R_offdiag_ > 0) {
        cuda_gmg_compute_restriction(R_diag_, col_R_diag_, row_R_diag_, nrows_,
                                     R_offdiag_, col_R_offdiag_, row_R_offdiag_,
                                     dev_res_, dev_res_ghost_, dev_tmp_,
                                     grid_dim_, block_dim_, stream_);
      } else {
        cuda_gmg_compute_restriction_no_offdiag(
            R_diag_, col_R_diag_, row_R_diag_, nrows_, dev_res_, dev_tmp_,
            grid_dim_, block_dim_, stream_);
      }
#endif
    }
  }

  void ComputeUpdatedSolutionOnDevice(void) {
    if (this->is_scheduled_to_this_process()) {
      assert(grid_transfer_on_device_);
#ifdef WITH_CUDA
      cudaSetDevice(0);

      if (nnz_P_offdiag_ > 0) {
        cuda_gmg_compute_updated_solution(
            P_diag_, col_P_diag_, row_P_diag_, nrows_, P_offdiag_,
            col_P_offdiag_, row_P_offdiag_, dev_sol_, dev_res_, dev_res_ghost_,
            grid_dim_, block_dim_, stream_);
      } else {
        cuda_gmg_compute_updated_solution_no_offdiag(
            P_diag_, col_P_diag_, row_P_diag_, nrows_, dev_sol_, dev_res_,
            grid_dim_, block_dim_, stream_);
      }
#endif
    }
  }

  void CopySolToHost(void) {
    if (this->is_scheduled_to_this_process()) {
      assert(grid_transfer_on_device_);
#ifdef WITH_CUDA
      cudaSetDevice(0);
      cudaMemcpyAsync(host_sol_, dev_sol_, nrows_ * sizeof(DataType),
                      cudaMemcpyDeviceToHost, stream_);
      cudaStreamSynchronize(stream_);
      memcpy(lvec_sol_->buffer, host_sol_, nrows_ * sizeof(DataType));
#endif
    }
  }

  void CopyResToHost(void) {
    if (this->is_scheduled_to_this_process()) {
      assert(grid_transfer_on_device_);
#ifdef WITH_CUDA
      cudaSetDevice(0);
      cudaMemcpyAsync(host_res_, dev_res_, nrows_ * sizeof(DataType),
                      cudaMemcpyDeviceToHost, stream_);
      cudaStreamSynchronize(stream_);
      memcpy(lvec_res_->buffer, host_res_, nrows_ * sizeof(DataType));
#endif
    }
  }

  void CopyResToTmp(void) {
    if (this->is_scheduled_to_this_process()) {
      assert(grid_transfer_on_device_);
#ifdef WITH_CUDA
      cudaSetDevice(0);
      cudaMemcpyAsync(dev_tmp_, dev_res_, nrows_ * sizeof(DataType),
                      cudaMemcpyDeviceToDevice, stream_);
#endif
    }
  }

  //   void prepare_interpolation_on_device(void);

  void clear(void) {
#ifdef WITH_CUDA
    if (A_diag_)
      cudaFree(A_diag_);
    if (A_offdiag_)
      cudaFree(A_offdiag_);
    if (R_diag_)
      cudaFree(R_diag_);
    if (R_offdiag_)
      cudaFree(R_offdiag_);
    if (P_diag_)
      cudaFree(P_diag_);
    if (P_offdiag_)
      cudaFree(P_offdiag_);
    if (async_iter_gpu_ || dev_sol_)
      cudaFree(dev_sol_);
    if (async_iter_gpu_ || host_sol_)
      cudaFreeHost(host_sol_);
    if (async_iter_gpu_ || dev_sol_ghost_)
      cudaFree(dev_sol_ghost_);
    if (async_iter_gpu_ || host_sol_ghost_)
      cudaFreeHost(host_sol_ghost_);
    if (async_iter_gpu_ || dev_rhs_)
      cudaFree(dev_rhs_);
    if (async_iter_gpu_ || host_rhs_)
      cudaFreeHost(host_rhs_);
    if (dev_rhs_ghost_)
      cudaFree(dev_rhs_ghost_);
    if (host_rhs_ghost_)
      cudaFreeHost(host_rhs_ghost_);
    if (dev_res_)
      cudaFree(dev_res_);
    if (host_res_)
      cudaFreeHost(host_res_);
    if (dev_res_ghost_)
      cudaFree(dev_res_ghost_);
    if (host_res_ghost_)
      cudaFreeHost(host_res_ghost_);
    if (dev_tmp_)
      cudaFree(dev_tmp_);
    if (row_A_diag_)
      cudaFree(row_A_diag_);
    if (row_A_offdiag_)
      cudaFree(row_A_offdiag_);
    if (row_R_diag_)
      cudaFree(row_R_diag_);
    if (row_R_offdiag_)
      cudaFree(row_R_offdiag_);
    if (row_P_diag_)
      cudaFree(row_P_diag_);
    if (row_P_offdiag_)
      cudaFree(row_P_offdiag_);
    if (col_A_diag_)
      cudaFree(col_A_diag_);
    if (col_A_offdiag_)
      cudaFree(col_A_offdiag_);
    if (col_R_diag_)
      cudaFree(col_R_diag_);
    if (col_R_offdiag_)
      cudaFree(col_R_offdiag_);
    if (col_P_diag_)
      cudaFree(col_P_diag_);
    if (col_P_offdiag_)
      cudaFree(col_P_offdiag_);
#endif
  }

private:
  VectorPtr res_;
  bool grid_transfer_on_device_;

#ifdef WITH_CUDA
  DataType *A_diag_, *A_offdiag_, *R_diag_, *R_offdiag_, *P_diag_, *P_offdiag_;
  DataType *dev_sol_, *dev_sol_ghost_, *dev_rhs_, *dev_rhs_ghost_, *dev_res_,
      *dev_res_ghost_, *dev_tmp_;
  DataType *host_sol_, *host_sol_ghost_, *host_rhs_, *host_rhs_ghost_,
      *host_res_, *host_res_ghost_;
  int *row_A_diag_, *row_A_offdiag_, *col_A_diag_, *col_A_offdiag_;
  int *row_R_diag_, *row_R_offdiag_, *col_R_diag_, *col_R_offdiag_;
  int *row_P_diag_, *row_P_offdiag_, *col_P_diag_, *col_P_offdiag_;
  int nnz_A_diag_, nnz_A_offdiag_, nnz_R_diag_, nnz_R_offdiag_, nnz_P_diag_,
      nnz_P_offdiag_;
  int nrows_;
  int grid_dim_, block_dim_;

  CPU_lVector< DataType > *lvec_sol_, *lvec_sol_ghost_, *lvec_rhs_,
      *lvec_rhs_ghost_, *lvec_res_, *lvec_res_ghost_;
  AsynchronousIterationGPU< LAD > *async_iter_gpu_;
  cudaStream_t stream_;
#endif
};

// function to be registered as signal handler for SIGUSR1, does nothing;
// default signal handler might terminate the process, but it shall just wake up

void do_nothing(int unused);

template < class LAD, int DIM >
BasicLevel< LAD, DIM >::BasicLevel(MPI_Comm global_comm, MPI_Comm partial_comm,
                              const MeshPtr &local_mesh, int master_rank,
                              const std::vector< int > &fe_degrees,
                              StandardGlobalAssembler< DataType, DIM > *global_asm,
                              const Settings &settings)
    : global_comm_(global_comm), partial_comm_(partial_comm),
      global_asm_(global_asm), fe_degrees_(fe_degrees),
      master_rank_(master_rank), settings_(settings), global_size_(0),
      global_rank_(-1), partial_rank_(-1) {
  if (global_comm_ != MPI_COMM_NULL) {
    assert(partial_comm_ != MPI_COMM_NULL);

    MPI_Comm_rank(global_comm_, &global_rank_);
    MPI_Comm_rank(partial_comm_, &partial_rank_);
    MPI_Comm_size(global_comm_, &global_size_);

    LOG_INFO("BasicLevel",
             "Creating level with " << global_size_ << " processes.");

    SharedVertexTable shared_verts;
    mesh_ = compute_ghost_cells(*local_mesh, global_comm_, shared_verts);

    init_la();

    // share pid of all procs of the partial comm
    int partial_size = 0;
    MPI_Comm_size(partial_comm_, &partial_size);
    assert(partial_size > 0);
    int my_pid = static_cast< int >(getpid());
    pids_ = std::vector< int >(partial_size, 0);
    MPI_Allgather(&my_pid, 1, MPI_INT, util::raw_array(pids_), 1, MPI_INT,
                  partial_comm_);
  }

  // register do_nothing function for signal SIGUSR1
  signal(SIGUSR1, do_nothing);
}

template < class LAD, int DIM >
void BasicLevel< LAD, DIM >::assemble_system(
    typename GlobalAssemblerType::VectorAssemblyFunction v_asm,
    typename GlobalAssemblerType::MatrixAssemblyFunction m_asm) {
  if (is_scheduled_to_this_process()) {
    assemble_matrix(m_asm);

    global_asm_->assemble_vector(*space_, v_asm, *rhs_);
  }
}

template < class LAD, int DIM >
void BasicLevel< LAD, DIM >::assemble_matrix(
    typename GlobalAssemblerType::MatrixAssemblyFunction m_asm) {
  if (is_scheduled_to_this_process()) {
    global_asm_->assemble_matrix(*space_, m_asm, *matrix_);
    matrix_->Compress();
  }
}

template < class LAD, int DIM >
void BasicLevel< LAD, DIM >::visualize_vector(const VectorPtr v,
                                         const std::string &filename) const {
  if (is_scheduled_to_this_process()) {
    v->Update();

    int max_fe_degree =
        *std::max_element(fe_degrees_.begin(), fe_degrees_.end());

    int num_intervals = max_fe_degree;
    CellVisualization< DataType, DIM > visu(*space_, num_intervals);

    std::vector< DataType > remote_index(mesh_->num_entities(mesh_->tdim()), 0);
    std::vector< DataType > sub_domain(mesh_->num_entities(mesh_->tdim()), 0);
    std::vector< DataType > material_number(mesh_->num_entities(mesh_->tdim()),
                                            0);

    for (mesh::EntityIterator it = mesh_->begin(mesh_->tdim());
         it != mesh_->end(mesh_->tdim()); ++it) {
      int temp1, temp2;
      //       mesh_->get_attribute_value("_remote_index_", mesh_->tdim(),
      //                                   it->index(),
      //                                   &temp1);
      mesh_->get_attribute_value("_sub_domain_", mesh_->tdim(), it->index(),
                                 &temp2);
      remote_index.at(it->index()) = temp1;
      sub_domain.at(it->index()) = temp2;
      material_number.at(it->index()) =
          mesh_->get_material_number(mesh_->tdim(), it->index());
    }
	std::vector<size_t> visu_vars;
	std::vector<std::string> names;
    for (int var = 0; var < space_->nb_var(); ++var) {
      std::string identifier = "var";
      identifier += boost::lexical_cast< std::string >(var);
	  visu_vars.push_back(var);	
	  names.push_back(identifier);
      
    }

    visu.visualize(FeEvalCell< DataType, DIM >(*space_, *v, visu_vars), names);



    //     visu.visualize_cell_data(remote_index, "_remote_index_");
    visu.visualize_cell_data(sub_domain, "_sub_domain_");
    visu.visualize_cell_data(material_number, "Material Id");
    

    VTKWriter< DataType, DIM> vtk_writer (visu, global_comm_, master_rank_);
    vtk_writer.write(filename);

    
  }
}

template < class LAD, int DIM >
typename BasicLevel< LAD, DIM >::VectorPtr
BasicLevel< LAD, DIM >::create_vector(void) const {
  if (is_scheduled_to_this_process()) {
    VectorPtr vec = VectorPtr(new VectorType());
    vec->Init(global_comm_, *couplings_, settings_.platform,
              settings_.implementation);
    vec->InitStructure();
    vec->Zeros();
    return vec;
  } else
    return VectorPtr();
}

template < class LAD, int DIM >
typename BasicLevel< LAD, DIM >::MatrixPtr
BasicLevel< LAD, DIM >::create_matrix(void) const {
  if (is_scheduled_to_this_process()) {
    MatrixPtr mat = MatrixPtr(new MatrixType());
    mat->Init(global_comm_, *couplings_, settings_.platform,
              settings_.implementation, settings_.matrix_format);

    SparsityStructure sparsity;
    global_asm_->compute_sparsity_structure(*space_, sparsity);

    mat->InitStructure(
        vec2ptr(sparsity.diagonal_rows), vec2ptr(sparsity.diagonal_cols),
        sparsity.diagonal_rows.size(), vec2ptr(sparsity.off_diagonal_rows),
        vec2ptr(sparsity.off_diagonal_cols), sparsity.off_diagonal_rows.size());
    mat->Zeros();
    return mat;
  } else
    return MatrixPtr();
}

template < class LAD, int DIM >
typename BasicLevel< LAD, DIM >::MatrixPtr
BasicLevel< LAD, DIM >::create_matrix(const LaCouplings &cp,
                                 const SparsityStructure &sp) const {
  if (is_scheduled_to_this_process()) {
    MatrixPtr mat = MatrixPtr(new MatrixType());
    mat->Init(global_comm_, cp, settings_.platform, settings_.implementation,
              settings_.matrix_format);
    mat->InitStructure(vec2ptr(sp.diagonal_rows), vec2ptr(sp.diagonal_cols),
                       sp.diagonal_rows.size(), vec2ptr(sp.off_diagonal_rows),
                       vec2ptr(sp.off_diagonal_cols),
                       sp.off_diagonal_rows.size());
    mat->Zeros();
    return mat;
  } else
    return MatrixPtr();
}

template < class LAD, int DIM > 
void BasicLevel< LAD, DIM >::create_restriction_matrix(void) {
  restriction_matrix_ = create_matrix();
  if (!tmp_)
    tmp_ = create_vector();
}

template < class LAD, int DIM >
void BasicLevel< LAD, DIM >::create_interpolation_matrix(void) {
  interpolation_matrix_ = create_matrix();
  if (!tmp_)
    tmp_ = create_vector();
}

template < class LAD, int DIM > 
void BasicLevel< LAD, DIM >::init_la() {
  space_ = SpacePtr(new VectorSpace< DataType, DIM >(global_comm_));
  couplings_ = CouplingsPtr(new LaCouplings());
 
  //TODO: fe_ansatz and is_cg hwo to choose??
  std::vector< bool > is_cg(space_->nb_var(),true);
  std::vector < FEType > fe_ansatz(space_->nb_var(), FEType::FE_TYPE_LAGRANGE);
  space_->Init(*mesh_, fe_ansatz, is_cg, fe_degrees_, hiflow::doffem::HIFLOW_CLASSIC);


  couplings_->Init(global_comm_);
  std::vector< int > global_offsets;
  std::vector< int > ghost_dofs;
  std::vector< int > ghost_offsets;
  space_->get_la_couplings(global_offsets, ghost_dofs, ghost_offsets);

  couplings_->InitializeCouplings(global_offsets, ghost_dofs, ghost_offsets);

  SparsityStructure sparsity;
  global_asm_->compute_sparsity_structure(*space_, sparsity);

  matrix_ = MatrixPtr(new MatrixType());

  matrix_->Init(global_comm_, *couplings_, settings_.platform,
                settings_.implementation, settings_.matrix_format);

  matrix_->InitStructure(
      vec2ptr(sparsity.diagonal_rows), vec2ptr(sparsity.diagonal_cols),
      sparsity.diagonal_rows.size(), vec2ptr(sparsity.off_diagonal_rows),
      vec2ptr(sparsity.off_diagonal_cols), sparsity.off_diagonal_rows.size());

  matrix_->Zeros();

  rhs_ = create_vector();
  sol_ = create_vector();
}

// template<class LAD, class ConnectionType>
// void GMGLevel<LAD, ConnectionType>::prepare_interpolation_on_device(void)
// {
//   LOG_ERROR("GMGLevel::prepare_interpolation_on_device() not yet
//   implemented!"); exit(-1);
// }

} // namespace gmg
} // namespace la
} // namespace hiflow

#endif
