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

#ifndef HIFLOW_ADAPTIVITY_FE_INTERPOLATION
#define HIFLOW_ADAPTIVITY_FE_INTERPOLATION

/// \author Philipp Gerstner

#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <sstream>
#include <utility>

#include <boost/function.hpp>
#include <mpi.h>
#include "common/log.h"
#include "common/vector_algebra.h"
#include "linear_algebra/la_descriptor.h"
#include "mesh/types.h"

enum FEInterpolationType { FE_INTER_LAGRANGE = 1, FE_INTER_INNERPROD = 2 };

// TODO: check whether is_cg = false works
namespace hiflow {

namespace la {
template <class LAD, class PreLAD> class LinearSolver;
template <class LAD> class CG;
template <class DataType> class SeqDenseMatrix;
class LaCouplings;
}

namespace mesh {
class Mesh;
class Entity;
}

namespace doffem {

}

template <class DataType, int DIM> class Element;
template <class DataType, int DIM> class StandardGlobalAssembler;
template <class DataType, int DIM> class VectorSpace;

///
/// \class FEInterpolation fe_interpolation.h
/// \brief abstract base class for interpolating FE functions from one space to
/// another one

template < class LAD, int DIM > class FEInterpolation {
  typedef typename LAD::DataType DataType;
  typedef typename LAD::VectorType VectorType;
  typedef Vec<DIM, DataType> Coord;

public:
  FEInterpolation();

  virtual ~FEInterpolation() { this->clear(); }

  /// \brief Initialize fe interpolation, this function involves the follwoing
  /// steps: <br> check whether input space is suitable for patching <br> build
  /// dof interpoaltion map from interpolating space to input space
  /// @param[in] from_space space to be interpolated from
  /// @param[in] to_space space to be interpolated to
  virtual void init(VectorSpace< DataType, DIM > *from_space,
                    VectorSpace< DataType, DIM > *to_space,
                    const VectorType &from_vector, bool reuse = true);

  /// \brief apply fe interpolation
  /// @param[in] from_vector vector to be interpolated
  /// @param[out] to_vector interpolated vector

  virtual void interpolate(const VectorType &from_vector,
                           VectorType &to_vector){};

  /// \brief set tolerance for identifying dofs by physical coordinates
  /// @param[in] eps tolerance

  virtual void set_dof_tolerance(DataType eps) 
  { 
    this->eps_ = eps; 
  }

  /// \brief get space to be interpolated from
  /// @return reference to space

  virtual VectorSpace< DataType, DIM > *get_from_space() {
    return this->from_space_;
  }

  /// \brief get space to be interpolated to
  /// @return reference to space

  virtual VectorSpace< DataType, DIM > *get_to_space() { return this->to_space_; }

  /// \brief clear all data structs
  virtual void clear();

  virtual bool is_initialized() const 
  { 
    return this->initialized_; 
  }

protected:
  /// \brief checks if input space is suitable for patch interpolation


  virtual bool check_space() const { return true;};

  /// pointers to input and output space
  VectorSpace< DataType, DIM > *from_space_;
  VectorSpace< DataType, DIM > *to_space_;

  /// pointers to input and output meshes
  mesh::MeshPtr from_mesh_;
  mesh::MeshPtr to_mesh_;

  /// FE degrees of interpolating space
  std::vector< int > degrees_;

  /// CG / DG flags for interpolating space
  std::vector< bool > is_cg_;

  /// number of variables in input space
  int nb_fe_;

  /// tolerance for dof identification
  DataType eps_;

  /// flag indicating whether patchinterpolation is initialized for given space
  bool initialized_;

  /// MPI rank
  int rank_;

  /// topological dimension of input space
  int tdim_;

  /// geometrical dimension of input space
  int gdim_;

  /// reuse or clean initialized interpolator
  bool reuse_;
};

///
/// \class FEInterpolationLagrange fe_interpolation.h
/// \brief class for interpolating Lagrange FE functions from one space to
/// another one \br interpolation works in the following way: given u_f (x) =
/// \sum_i a_i \phi_i(x) \in V_f (from_space) \br compute u_t (x) = \sum_j b_j
/// \psi_j(x) \in V_t (to_space) such that \bre u_t(x_j) = u_f(x_j) for all dof
/// positions x_j (t-> to_space) with \psi_k(x_j) = delta_{kj}

template < class LAD, int DIM >
class FEInterpolationLagrange : virtual public FEInterpolation< LAD, DIM > {
  typedef typename LAD::DataType DataType;
  typedef typename LAD::VectorType VectorType;
  typedef Vec<DIM, DataType> Coord;

public:
  FEInterpolationLagrange();

  ~FEInterpolationLagrange() { this->clear(); }

  /// \brief Initialize fe interpolation, this function involves the follwoing
  /// steps: <br> check whether input space is suitable for patching <br> build
  /// dof interpoaltion map from interpolating space to input space
  /// @param[in] from_space space to be interpolated from
  /// @param[in] to_space space to be interpolated to
  void init(VectorSpace< DataType, DIM > *from_space,
            VectorSpace< DataType, DIM > *to_space, const VectorType &from_vector,
            bool reuse = true);

  /// \brief apply fe interpolation
  /// @param[in] from_vector vector to be interpolated
  /// @param[out] to_vector interpolated vector
  void interpolate(const VectorType &from_vector, VectorType &to_vector);

  /// \brief clear all data structs
  void clear();

protected:
  /// \brief checks if input space is suitable for patch interpolation
  bool check_space() const;

  /// \brief build dof interpolation mapping
  void create_dof_mapping();

  /// \brief dof_weights(j) = \{ (i_1, w_1), ... (i_n(j), w_n(j) \} such that
  /// \br b_j = \sum_{l} a_dof_weight(j).i_l * dof_weights(j).w_l
  std::vector< std::map< int, std::map< int, DataType > > > dof_weights_;
};

#if 0 // TODO (Philipp G.): debugging required
///
/// \class FEInterpolationInnerProduct fe_interpolation.h
/// \brief class for interpolating Lagrange FE functions from one space to
/// another one \br interpolation works in the following way: given u_f (x) =
/// \sum_i a_i \phi_i(x) \in V_f (from_space) \br compute u_t (x) = \sum_j b_j
/// \psi_j(x) \in V_t (to_space) such that \bre
///
template < class LAD, int DIM > class InnerProductAssembler;

template < class LAD, int DIM >
class FEInterpolationInnerProduct : virtual public FEInterpolation< LAD, DIM > 
{
  typedef typename LAD::DataType DataType;
  typedef typename LAD::VectorType VectorType;
  typedef typename LAD::MatrixType OperatorType;

public:
  FEInterpolationInnerProduct() : FEInterpolation< LAD, DIM >() 
  {
    this->comm_ = MPI_COMM_NULL;
  }

  ~FEInterpolationInnerProduct() 
  { 
    this->clear(); 
  }

  void set_inner_product_assembler(
      InnerProductAssembler< LAD, DIM > *local_asm,
      const std::vector< std::vector< bool > > &coupling_vars) {
    assert(local_asm != nullptr);
    assert(coupling_vars.size() > 0);

    this->local_asm_ = local_asm;
    this->to_coupling_vars_ = coupling_vars;
  }

  void set_dirichlet_bc(std::vector< int > dofs, std::vector< DataType > values) 
  {
    this->to_dirichlet_dofs_ = dofs;
    this->to_dirichlet_values_ = values;
  }

  void set_solver(la::LinearSolver< LAD > *solver) { this->user_solver_ = solver; }

  /// \brief Initialize fe interpolation, this function involves the follwoing
  /// steps: <br> check whether input space is suitable for patching <br> build
  /// dof interpoaltion map from interpolating space to input space
  /// @param[in] from_space space to be interpolated from
  /// @param[in] to_space space to be interpolated to
  void init(VectorSpace< DataType, DIM > *from_space,
            VectorSpace< DataType, DIM > *to_space, const VectorType &from_vector,
            bool reuse = true) {
    NOT_YET_IMPLEMENTED;
  }

  /// \brief apply fe interpolation
  /// @param[in] from_vector vector to be interpolated
  /// @param[out] to_vector interpolated vector
  void interpolate(const VectorType &from_vector, VectorType &to_vector) 
  {
    NOT_YET_IMPLEMENTED;
  }

  /// \brief clear all data structs
  void clear() 
  {
    FEInterpolation< LAD, DIM >::clear();
    // this->local_asm_ = nullptr;
    // this->to_coupling_vars_.clear();
    this->mass_matrix_.Clear();
    this->rhs_.Clear();
    this->cg_.Clear();
    // this->amg_.Clear();
    // this->ilupp_.Clear();
    // this->to_dirichlet_dofs_.clear();
    // this->to_dirichlet_values_.clear();
    // this->user_solver_ = nullptr;
    this->to_dof_func_.clear();
    // this->from_vector_ = nullptr;
  }

  VectorType *get_rhs() { return &this->rhs_; }

protected:
  void setup_asm() { NOT_YET_IMPLEMENTED; }

  void setup_LA() { NOT_YET_IMPLEMENTED; }

  void setup_solver() { NOT_YET_IMPLEMENTED; }

  /// \brief checks if input space is suitable for patch interpolation
  bool check_space() const { NOT_YET_IMPLEMENTED; return true;}

  /// \brief build dof interpolation mapping
  void assemble_matrix() { NOT_YET_IMPLEMENTED; }
  void assemble_rhs() { NOT_YET_IMPLEMENTED; }

  OperatorType mass_matrix_;
  VectorType rhs_;
  LaCouplings couplings_;

  la::LinearSolver< LAD > *user_solver_;

  la::CG< LAD > cg_;
  // PreconditionerIlupp<LAD> ilupp_;
  // HypreBoomerAMG<LAD> amg_;

  StandardGlobalAssembler< DataType, DIM > global_asm_;
  InnerProductAssembler< LAD, DIM > *local_asm_;

  std::vector< std::vector< bool > > to_coupling_vars_;

  MPI_Comm comm_;

  std::vector< int > to_dirichlet_dofs_;
  std::vector< DataType > to_dirichlet_values_;

  std::vector< int > to_dof_func_;

  VectorType const *from_vector_;
};
#endif

template < class LAD, int DIM > 
class InnerProductAssembler {
  typedef typename LAD::DataType DataType;
  typedef typename LAD::VectorType VectorType;
  typedef la::SeqDenseMatrix< DataType > LocalMatrix;
  typedef std::vector< DataType > LocalVector;
  
public:
  enum SIDE
  {
    LEFT = 0,
    RIGHT = 1
  };
  
  InnerProductAssembler()
  : left_initialized_(false),
  right_initialized_(false),
  left_vec_(nullptr),
  right_vec_ (nullptr),
  left_space_(nullptr),
  right_space_ (nullptr)
  {}

  ~InnerProductAssembler() {}

  virtual void operator()(const Element< DataType, DIM > &element,
                          const Quadrature< DataType > &quadrature,
                          LocalMatrix &lm) = 0;

  virtual void operator()(const Element< DataType, DIM > &element,
                          const Quadrature< DataType > &quadrature,
                          LocalVector &lv) = 0;

  virtual void operator()(const Element< DataType > &element,
                          const Quadrature< DataType > &quadrature,
                          DataType &ls) = 0;
                          
  void set_vector(SIDE side, VectorType const *vector) 
  {
    assert(vector != nullptr);
    switch (side)
    {
      case LEFT:
        this->left_vec_ = vector;
        this->left_initialized_ = false;
        break;
      case RIGHT:
        this->right_vec_ = vector;
        this->right_initialized_ = false;
        break;
      default:
        assert(0);
        break;
    }
  }

  void set_space(SIDE side, VectorSpace< DataType, DIM > const *space) 
  {
    assert(space != nullptr);
    switch (side)
    {
      case LEFT:
        this->left_space_ = space;
        this->left_initialized_ = false;
        break;
      case RIGHT:
        this->right_space_ = space;
        this->right_initialized_ = false;
        break;
      default:
        assert(0);
        break;
    }
  }

protected:
  VectorType const &vector_left() const 
  { 
    return *this->left_vec_; 
  }

  VectorType const &vector_right() const 
  { 
    return *this->right_vec_; 
  }
  
  VectorType const *left_vec_;
  VectorType const *right_vec_;
  
  VectorSpace< DataType, DIM > const *left_space_;
  VectorSpace< DataType, DIM > const *right_space_;
  
  bool left_initialized_;
  bool right_initialized_;
};



} // namespace hiflow
#endif
