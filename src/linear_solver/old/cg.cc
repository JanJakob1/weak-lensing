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

/// @author Chandramowli Subramanian

#include <cassert>
#include <cmath>

#include "linear_solver/cg.h"
#include "common/log.h"
#include "linear_algebra/pce_matrix.h"
#include "linear_algebra/block_matrix.h"
#include "linear_algebra/la_descriptor.h"
#include <iomanip>

namespace hiflow {
namespace la {

/// standard constructor

template < class LAD > CG< LAD >::CG() : LinearSolver< LAD >() {
  this->SetMethod("NoPreconditioning");
  if (this->print_level_ > 2) {
    LOG_INFO("Linear solver", "CG");
    LOG_INFO("Preconditioning", this->Method());
  }
  this->name_ = "CG";
}

/// destructor

template < class LAD > CG< LAD >::~CG() {}

/// Sets parameters of the solution process
/// @param method "NoPreconditioning" or "Preconditioning" -- whether to use
/// preconditioning or not.

template < class LAD >
void CG< LAD >::InitParameter(std::string const &method) {
  // chose method_
  this->precond_method_ = method;
  assert((this->Method() == "NoPreconditioning") ||
         (this->Method() == "Preconditioning") ||
	 (this->Method() == "RightPreconditioning") ||
	 (this->Method() == "LeftPreconditioning"));
}

/// Solves the linear system.
/// @param [in] b right hand side vector
/// @param [in,out] x start and solution vector

template < class LAD >
LinearSolverState CG< LAD >::SolveImpl(const VectorType &b, VectorType *x) {
  assert(x->is_initialized());
  assert(b.is_initialized());

  LinearSolverState state;

  if (this->Method() == "NoPreconditioning") {
    state = this->SolveNoPrecond(b, x);
  } else if (this->Method() == "Preconditioning") {
    state = this->SolvePrecond(b, x);
  } else if (this->Method() == "RightPreconditioning") {
    state = this->SolvePrecond(b, x);
  } else if (this->Method() == "LeftPreconditioning") {
    state = this->SolvePrecond(b, x);
  } else {
    state = kSolverError;
  }

  return state;
}

/// Solve without preconditioning.

template < class LAD >
LinearSolverState CG< LAD >::SolveNoPrecond(const VectorType &b,
                                            VectorType *x) {
  assert(this->Method() == "NoPreconditioning");
  assert(this->op_ != nullptr);

  if (this->print_level_ > 2) {
    LOG_INFO(this->name_, "solve without preconditioning");
  }

  IterateControl::State conv = IterateControl::kIterate;

  // needed vectors
  VectorType r, p, Ap;
  r.CloneFromWithoutContent(b);
  p.CloneFromWithoutContent(b);
  Ap.CloneFromWithoutContent(b);

  // needed values
  DataType alpha, beta;

  // initialization step
  this->iter_ = 0;

  this->op_->VectorMult(*x, &r);

  r.ScaleAdd(b, static_cast< DataType >(-1.));

  DataType ressquared = r.Dot(r);
  this->res_ = std::sqrt(ressquared);
  this->res_init_ = this->res_;
  this->res_rel_ = 1.;
  conv = this->control().Check(this->iter_, this->res_);

  if (this->print_level_ > 1) {
    LOG_INFO(this->name_, "initial res norm   =  " << this->res_);
  }

  p.CopyFrom(r);
  beta = ressquared;

  // main loop
  while (conv == IterateControl::kIterate) {
    ++(this->iter_);
    this->op_->VectorMult(p, &Ap);
    alpha = beta / (Ap.Dot(p));
    x->Axpy(p, alpha);
    r.Axpy(Ap, -alpha);

    ressquared = r.Dot(r);
    this->res_ = sqrt(ressquared);
    this->res_rel_ = this->res_ / this->res_init_;
    if (this->print_level_ > 2) {
      LOG_INFO(this->name_,
               "residual (iteration " << this->iter_ << "): " << this->res_);
    }

    conv = this->control().Check(this->iter_, this->res_);
    if (conv != IterateControl::kIterate) {
      break;
    }

    beta = ressquared / beta;
    p.ScaleAdd(r, beta);
    beta = ressquared;
  }

  if (this->print_level_ > 1) {
    LOG_INFO(this->name_, "final iterations   = " << this->iter_);
    LOG_INFO(this->name_, "final abs res norm = " << this->res_);
    LOG_INFO(this->name_, "final rel res norm = " << this->res_rel_)
  } 

  if (conv == IterateControl::kFailureDivergenceTol ||
      conv == IterateControl::kFailureMaxitsExceeded) {
    return kSolverExceeded;
  }
  return kSolverSuccess;
}

/// Solve with preconditioning.

template < class LAD >
LinearSolverState CG< LAD >::SolvePrecond(const VectorType &b, VectorType *x) {
  assert(this->Method() != "NoPreconditioning");
  assert(this->op_ != nullptr);
  assert(this->precond_ != nullptr);

  if (this->print_level_ > 2) {
    LOG_INFO(this->name_, "solve with right preconditioning");
  }

  IterateControl::State conv = IterateControl::kIterate;

  // needed vectors
  VectorType r, p, z, Ap;
  r.CloneFromWithoutContent(b);
  p.CloneFromWithoutContent(b);
  z.CloneFromWithoutContent(b);
  Ap.CloneFromWithoutContent(b);

  // needed values
  DataType alpha, beta, gamma, ressquared;

  // initialization step
  this->iter_ = 0;
  std::cout << "test"<< std::endl;
  this->op_->VectorMult(*x, &r);
  std::cout << "test2"<< std::endl;
  std::vector<int> id;
  std::vector<DataType> val;
  r.GetAllDofsAndValues(id, val);
  std::cout << "Ax size: " << val.size() << std::endl;
  /*for (int i = 0; i < val.size(); ++i){
    std::cout << " Ax[" << i << "] " << val[i] << std::endl;
  }*/
  r.ScaleAdd(b, static_cast< DataType >(-1.));
  r.GetAllDofsAndValues(id, val);
  for (int i = 0; i < val.size(); ++i){
    std::cout << "Ax-b[" << i << "] " << val[i] << std::endl;
  }
  ressquared = r.Dot(r);
  this->res_init_ = this->res_ = sqrt(ressquared);
  this->res_rel_ = 1.;
  conv = this->control().Check(this->iter_, this->res_);
  z.Zeros();
  this->ApplyPreconditioner(r, &z);
  p.CopyFrom(z);

  if (this->print_level_ > 1) {
    LOG_INFO(this->name_, "initial res norm   = " << this->res_);
  }

  beta = r.Dot(z);

  // main loop
  while (conv == IterateControl::kIterate) {
    ++(this->iter_);
    this->op_->VectorMult(p, &Ap);
    alpha = beta / (Ap.Dot(p));
    x->Axpy(p, alpha);
    r.Axpy(Ap, -alpha);

    ressquared = r.Dot(r);
    this->res_ = sqrt(ressquared);
    this->res_rel_ = this->res_ / this->res_init_;
    if (this->print_level_ > 2) {
      LOG_INFO(this->name_,
               "residual (iteration " << this->iter_ << "): " << this->res_);
    }

    conv = this->control().Check(this->iter_, this->res_);
    if (conv != IterateControl::kIterate) {
      break;
    }

    z.Zeros();
    this->ApplyPreconditioner(r, &z);
    gamma = r.Dot(z);
    beta = gamma / beta;
    p.ScaleAdd(z, beta);
    beta = gamma;
  }

  if (this->print_level_ > 1) {
    LOG_INFO(this->name_, "final iterations   = " << this->iter_);
    LOG_INFO(this->name_, "final abs res norm = " << this->res_);
    LOG_INFO(this->name_, "final rel res norm = " << this->res_rel_);
  } 

  if (conv == IterateControl::kFailureDivergenceTol ||
      conv == IterateControl::kFailureMaxitsExceeded) {
    return kSolverExceeded;
  }
  return kSolverSuccess;
}

/// template instantiation
template class CG< LADescriptorCoupledD >;
template class CG< LADescriptorCoupledS >;
template class CG< LADescriptorPCE < LADescriptorCoupledD > >;
template class CG< LADescriptorPCE < LADescriptorCoupledS > >;
template class CG< LADescriptorBlock< LADescriptorCoupledD > >;
template class CG< LADescriptorBlock< LADescriptorCoupledS > >;
template class CG< LADescriptorPolynomialChaosD >;

#ifdef WITH_HYPRE
template class CG< LADescriptorHypreD >;
template class CG< LADescriptorPCE < LADescriptorHypreD > >;
template class CG< LADescriptorBlock< LADescriptorHypreD > >;
#endif

} // namespace la
} // namespace hiflow
