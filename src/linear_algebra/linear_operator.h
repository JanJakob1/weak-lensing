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

/// @author Philipp Gerstner

#ifndef HIFLOW_LINEAR_ALGEBRA_LINEAR_OPERATOR_H
#define HIFLOW_LINEAR_ALGEBRA_LINEAR_OPERATOR_H

#include "common/log.h"
#include "linear_algebra/vector.h"
#include "space/vector_space.h"
#include <assert.h>
#include <cstddef>

namespace hiflow {
namespace la {

/// \brief Abstract base class for distributed linear operator implementations.

template < class DataType > class LinearOperator {
public:
  /// Standard Constructor

  LinearOperator(){};

  /// Destructor

  virtual ~LinearOperator() {}
                         
  virtual void VectorMult(Vector< DataType > &in,
                          Vector< DataType > *out) const = 0;

  virtual void VectorMultOffdiag(Vector< DataType > &in,
                                 Vector< DataType > *out)
  {
    NOT_YET_IMPLEMENTED;
  }
                          
  /// out = beta * out + alpha * this * in

  void VectorMultAdd(DataType alpha, Vector< DataType > &in,
                     DataType beta, Vector< DataType > *out) const {
  
    assert(out != nullptr);
    Vector< DataType > *tmp = out->Clone();
    this->VectorMult(in, tmp);
    out->Scale(beta);
    out->Axpy(*tmp, alpha);
  }

  virtual bool IsInitialized() const { return true; }
};

template < class LAD, int DIM, class Application> 
class MatrixFree : public virtual LinearOperator<typename LAD::DataType> {
public:
  typedef typename LAD::VectorType VectorType;
  typedef typename LAD::DataType DataType;
  
  MatrixFree () :ap_(nullptr), space_(nullptr) {};
  
  ~MatrixFree  () {};
  
  void set_application(Application * ap) {
    assert(ap!=nullptr);
    this->ap_ = ap;
  }
  
  void set_space(VectorSpace<DataType, DIM>* space) {
    assert(space!= nullptr);
    this->space_ = space;    
  }
  
   void VectorMult(Vector< DataType >& in, Vector< DataType > *out) const{
    out->Update(); 
    assert(this->ap_!=nullptr);
    assert(this->space_!=nullptr);
    in.Update();
    this->ap_->VectorMult(space_, &in, out);
    out->Update();
  }
  Application * ap_;
  VectorSpace<DataType, DIM>* space_;
};

template < class DataType >
class ZeroOperator : public LinearOperator< DataType > {
public:
  /// Standard Constructor

  ZeroOperator() { this->is_initialized_ = true; };
  /// Destructor

  virtual ~ZeroOperator() {}

  /// out = this * in

  virtual void VectorMult(Vector< DataType > &in,
                          Vector< DataType > *out) const {
    assert(out != nullptr);
    assert(this->IsInitialized());
    out->Zeros();
  }
};

/// \brief class for getting a chained operator out of two individual ones

template < class DataType >
class ChainLinearOperator : public LinearOperator< DataType > {
public:
  /// Standard Constructor

  ChainLinearOperator(){};
  /// Destructor

  virtual ~ChainLinearOperator() {}

  virtual void SetOperatorA(const LinearOperator< DataType > &op) {
    this->op_A_ = &op;
  }

  virtual void SetOperatorB(const LinearOperator< DataType > &op) {
    this->op_B_ = &op;
  }

  /// out = B[ A [n] ]

  virtual void VectorMult(Vector< DataType > &in,
                          Vector< DataType > *out) const {
    assert(op_A_ != nullptr);
    assert(op_B_ != nullptr);
    assert(out != nullptr);
    assert(this->IsInitialized());

    Vector< DataType > *tmp = in.CloneWithoutContent();

    // tmp = A[in]
    this->op_A_->VectorMult(in, tmp);

    // out = b[tmp]
    this->op_B_->VectorMult(*tmp, out);
  }

  virtual bool IsInitialized() const {
    if (this->op_A_ != nullptr) {
      if (this->op_A_->IsInitialized()) {
        if (this->op_B_ != nullptr) {
          if (this->op_B_->IsInitialized()) {
            return true;
          } else {
            LOG_DEBUG(0, "op_B not initialized");
          }
        } else {
          LOG_DEBUG(0, "op_B = nullptr");
        }
      } else {
        LOG_DEBUG(0, "op_A not initialized");
      }
    } else {
      LOG_DEBUG(0, "op_A = nullptr");
    }
    return false;
  }

protected:
  LinearOperator< DataType > const *op_A_;
  LinearOperator< DataType > const *op_B_;
};

/// \brief class for getting a linear combination of two linear operators

template < class DataType >
class SumLinearOperator : public LinearOperator< DataType > {
public:
  /// Standard Constructor

  SumLinearOperator(){};
  /// Destructor

  virtual ~SumLinearOperator() {}

  virtual void SetOperatorA(const LinearOperator< DataType > &op,
                            const DataType scale) {
    assert(op != nullptr);
    this->op_A_ = &op;
    this->scale_A_ = scale;
  }

  virtual void SetOperatorB(const LinearOperator< DataType > &op,
                            const DataType scale) {
    assert(op != nullptr);
    this->op_B_ = &op;
    this->scale_B = scale;
  }

  /// out = scale_A  * A[in] + scale_B * B[in]

  virtual void VectorMult(Vector< DataType > &in,
                          Vector< DataType > *out) const {
    assert(op_A_ != nullptr);
    assert(op_B_ != nullptr);
    assert(out != nullptr);
    assert(this->IsInitialized());

    // out = A[in]
    this->op_A_->VectorMult(in, out);

    // out = scale_A * out + scale_B * B[in]
    this->op_B_->VectorMultAdd(scale_B_, in, scale_A_, out);
  }

  virtual bool IsInitialized() const {
    if (this->op_A_ != nullptr) {
      if (this->op_A_->IsInitialized()) {
        if (this->op_B_ != nullptr) {
          if (this->op_B_->IsInitialized()) {
            return true;
          } else {
            LOG_DEBUG(0, "op_B not initialized");
          }
        } else {
          LOG_DEBUG(0, "op_B = nullptr");
        }
      } else {
        LOG_DEBUG(0, "op_A not initialized");
      }
    } else {
      LOG_DEBUG(0, "op_A = nullptr");
    }
    return false;
  }

protected:
  LinearOperator< DataType > const *op_A_;
  LinearOperator< DataType > const *op_B_;

  DataType scale_A_;
  DataType scale_B_;
};

} // namespace la
} // namespace hiflow

#endif
