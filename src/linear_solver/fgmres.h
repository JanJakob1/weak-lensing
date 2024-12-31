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

/// @author Hendryk Bockelmann, Chandramowli Subramanian

#ifndef HIFLOW_LINEARSOLVER_FGMRES_H_
#define HIFLOW_LINEARSOLVER_FGMRES_H_

#include "linear_solver/gmres.h"
#include "linear_solver/linear_solver.h"
#include <iomanip>
#include <string>
#include <vector>

namespace hiflow {
namespace la {

template < class DataType > class SeqDenseMatrix;

/// @brief Flexible GMRES solver
///
/// Flexible GMRES solver for linear systems Ax=b.

template < class LAD, class PreLAD = LAD >
class FGMRES : public GMRES< LAD, PreLAD > {
public:
  typedef typename LAD::MatrixType OperatorType;
  typedef typename LAD::VectorType VectorType;
  typedef typename LAD::DataType DataType;

  FGMRES();
  virtual ~FGMRES();

  virtual void Clear() {
    LinearSolver< LAD, PreLAD >::Clear();
    this->name_ = "FGMRES";
  }

  /// Initialize pointers for Krylov Basis
  virtual void InitBasis(const VectorType &ref_vec, int iteration);

  /// Allocate Krylov basis vectors
  virtual void AllocateBasis(int basis_size);

  /// Set all basis vectors to zero
  virtual void SetBasisToZero();

  /// Deallocate Krylov basis vectors
  virtual void FreeBasis();

private:
  LinearSolverState SolveLeft(const VectorType &b, VectorType *x);
  LinearSolverState SolveRight(const VectorType &b, VectorType *x);

  std::vector< VectorType * > Z_;
};

template < class LAD, class PreLAD >
void setup_FGMRES_solver(FGMRES< LAD, PreLAD > &fgmres_solver,
                         const PropertyTree &params,
                         NonlinearProblem< LAD > *nonlin) 
{
  const int max_it = params["MaxIt"].get< int >(1000);
  const int max_size = params["KrylovSize"].get< int >(100);
  const double abs_tol = params["AbsTol"].get< double >(1e-12);
  const double rel_tol = params["RelTol"].get< double >(1e-10);
  const bool use_press_filter = params["UsePressureFilter"].get< bool >(false);

  fgmres_solver.InitControl(max_it, abs_tol, rel_tol, 1e6);
  if (params["UsePrecond"].get< bool >(true)) 
  {
    fgmres_solver.InitParameter(max_size, "RightPreconditioning");
  } 
  else 
  {
    fgmres_solver.InitParameter(max_size, "NoPreconditioning");
  }
  fgmres_solver.SetPrintLevel(params["PrintLevel"].get< int >(0));
  fgmres_solver.SetReuse(params["Reuse"].get< bool >(true));
  fgmres_solver.SetReuseBasis(params["ReuseBasis"].get< bool >(true));

  if (use_press_filter && nonlin != nullptr) 
  {
    fgmres_solver.SetupNonLinProblem(nonlin);
  }
  fgmres_solver.SetName(params["Name"].get< std::string >("FGMRES"));
}

template < class LAD >
void setup_FGMRES_solver(FGMRES< LAD > &fgmres_solver,
                         const PropertyTree &params,
                         NonlinearProblem< LAD > *nonlin) {
  setup_FGMRES_solver< LAD, LAD >(fgmres_solver, params, nonlin);
}

/// @brief GMRES creator class
/// @author Tobias Hahn

template < class LAD > class FGMREScreator : public LinearSolverCreator< LAD > {
public:
  LinearSolver< LAD > *params(const PropertyTree &c) {
    FGMRES< LAD > *newFGMRES = new FGMRES< LAD >();
    if (c.contains("Method") && c.contains("SizeBasis")) {
      newFGMRES->InitParameter(
          c["SizeBasis"].template get< int >(),
          c["Method"].template get< std::string >().c_str());
    }
    if (c.contains("MaxIterations") && c.contains("AbsTolerance") &&
        c.contains("RelTolerance") && c.contains("DivTolerance")) {
      newFGMRES->InitControl(c["MaxIterations"].template get< int >(),
                            c["AbsTolerance"].template get< double >(),
                            c["RelTolerance"].template get< double >(),
                            c["DivTolerance"].template get< double >());
    }
    return newFGMRES;
  }
};

} // namespace la
} // namespace hiflow

#endif // HIFLOW_LINEARSOLVER_FGMRES_H_
