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

/// @author Dimitar Lukarski, Nico Trost, Niels Wegh

#include "linear_solver/preconditioner_ai.h"
#include "linear_algebra/lmp/lpreconditioner.h"
#include "linear_algebra/lmp/lpreconditioner_ai.h"
#include "linear_solver/preconditioner.h"

#include <cassert>

#include "common/log.h"
#include "linear_algebra/la_descriptor.h"

namespace hiflow {
namespace la {

template < class LAD >
PreconditionerApproximateInverse< LAD >::PreconditionerApproximateInverse()
    : Preconditioner< LAD >() {
  this->local_ai_ = nullptr;
}

template < class LAD >
PreconditionerApproximateInverse< LAD >::~PreconditionerApproximateInverse() {
  this->Clear();
}

template < class LAD >
void PreconditionerApproximateInverse< LAD >::Init_FSAI(const int power) {
  this->Clear();
  lPreconditioner_ApproximateInverse_FSAI< typename LAD::DataType > *lp =
      new lPreconditioner_ApproximateInverse_FSAI< typename LAD::DataType >;

  lp->Init(); // default values
  lp->set_matrix_power(power);
  this->local_ai_ = lp;
}

template < class LAD >
void PreconditionerApproximateInverse< LAD >::SetupOperator(OperatorType &op) {
  this->op_ = &op;
  this->local_ai_->SetupOperator(op.diagonal());

  this->SetModifiedOperator(true);
  this->local_ai_->SetModifiedOperator(true);
}

template < class LAD > void PreconditionerApproximateInverse< LAD >::Clear() {
  if (this->local_ai_ != nullptr) {
    this->local_ai_->Clear();
    delete this->local_ai_;
  }
  this->local_ai_ = nullptr;
  Preconditioner< LAD >::Clear();
}

template < class LAD >
void PreconditionerApproximateInverse< LAD >::SetupVector(
    const VectorType &vec) {
  this->local_ai_->SetupVector(&vec.interior());
}

template < class LAD >
void PreconditionerApproximateInverse< LAD >::BuildImpl(VectorType const *b,
                                                        VectorType *x) {
  assert(this->local_ai_ != nullptr);

  this->local_ai_->Build();
  this->local_ai_->SetModifiedOperator(false);
  this->local_ai_->SetState(true);
}

template < class LAD >
LinearSolverState
PreconditionerApproximateInverse< LAD >::SolveImpl(const VectorType &b,
                                                   VectorType *x) {
  assert(this->op_ != nullptr);
  assert(this->local_ai_ != nullptr);

  this->local_ai_->ApplylPreconditioner(b.interior(), &(x->interior()));

  return kSolverSuccess;
}

template < class LAD >
void PreconditionerApproximateInverse< LAD >::Print(std::ostream &out) const {
  this->local_ai_->print(out);
}

/// template instantiation
template class PreconditionerApproximateInverse< LADescriptorCoupledD >;
template class PreconditionerApproximateInverse< LADescriptorCoupledS >;

} // namespace la
} // namespace hiflow
