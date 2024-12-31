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

#ifndef HIFLOW_NONLINEAR_NONLINEAR_SOLVER_CREATOR_H_
#define HIFLOW_NONLINEAR_NONLINEAR_SOLVER_CREATOR_H_

#include "common/property_tree.h"
#include "nonlinear/nonlinear_solver.h"

namespace hiflow {

/// @brief Creator base class for linear solvers in HiFlow.
/// @author Tobias Hahn

template < class LAD, int DIM > class NonlinearSolverCreator {
public:
  typedef typename LAD::MatrixType MatrixType;
  typedef typename LAD::VectorType VectorType;

  virtual NonlinearSolver< LAD, DIM > *
  params(VectorType *residual, MatrixType *matrix, const PropertyTree &c) = 0;
  virtual NonlinearSolver< LAD, DIM > *params(const PropertyTree &c) = 0;
};

} // namespace hiflow

#endif // HIFLOW_NONLINEAR_NONLINEAR_SOLVER_CREATOR_H_
