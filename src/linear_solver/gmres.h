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

#ifndef HIFLOW_LINEARSOLVER_GMRES_H_
#define HIFLOW_LINEARSOLVER_GMRES_H_

#include "common/timer.h"
#include "linear_algebra/seq_dense_matrix.h"
#include "linear_solver/linear_solver.h"
#include "linear_solver/linear_solver_creator.h"
#include <cmath>
#include <string>
#include <vector>

namespace hiflow {
namespace la {

/// @brief GMRES solver
///
/// GMRES solver for linear systems Ax=b with left, right or no preconditioning.
/// (not flexible!)

template < class LAD, class PreLAD = LAD >
class GMRES : public LinearSolver< LAD, PreLAD > {
public:
  typedef typename LAD::MatrixType OperatorType;
  typedef typename LAD::VectorType VectorType;
  typedef typename LAD::DataType DataType;

  GMRES();
  virtual ~GMRES();

  virtual void InitParameter(int size_basis, std::string method);

  virtual int size_basis() const { return this->size_basis_; }

  virtual void Clear() {
    LinearSolver< LAD, PreLAD >::Clear();
    this->name_ = "GMRES";
  }

  /// Set flag whether or not Krylov basis should be reused

  virtual void SetReuseBasis(bool flag) { this->reuse_basis_ = flag; }
  /// Sets the relative tolerance.
  /// Needed by Inexact Newton Methods
  /// @param reltol relative tolerance of residual to converge

  virtual void SetRelativeTolerance(double reltol) {
    int maxits = this->control_.maxits();
    double atol = this->control_.absolute_tol();
    double dtol = this->control_.divergence_tol();
    this->control_.Init(maxits, atol, reltol, dtol);
  }

  /// Initialize pointers for Krylov Basis
  virtual void InitBasis(const VectorType &ref_vec, int iteration);

  /// Allocate Krylov basis vectors
  virtual void AllocateBasis(int basis_size);

  /// Set all basis vectors to zero
  virtual void SetBasisToZero();

  /// Deallocate Krylov basis vectors
  virtual void FreeBasis();

protected:
  /// Applies Givens rotation.
  /// @param cs cos(phi)
  /// @param sn sin(phi)
  /// @param dx first coordinate
  /// @param dy second coordinate

  inline virtual void ApplyPlaneRotation(const DataType &cs, const DataType &sn,
                                         DataType *dx, DataType *dy) const {
    const DataType temp = cs * (*dx) + sn * (*dy);
    *dy = -sn * (*dx) + cs * (*dy);
    *dx = temp;
  }

  /// Generates Givens rotation.
  /// @param dx first coordinate
  /// @param dy second coordinate
  /// @param cs cos(phi)
  /// @param sn sin(phi)

  inline virtual void GeneratePlaneRotation(const DataType &dx,
                                            const DataType &dy, DataType *cs,
                                            DataType *sn) const {
    const DataType beta = std::sqrt(dx * dx + dy * dy);
    *cs = dx / beta;
    *sn = dy / beta;
  }

  virtual LinearSolverState SolveImpl(const VectorType &b, VectorType *x);

  virtual LinearSolverState SolveNoPrecond(const VectorType &b, VectorType *x);
  virtual LinearSolverState SolveLeft(const VectorType &b, VectorType *x);
  virtual LinearSolverState SolveRight(const VectorType &b, VectorType *x);

  virtual void UpdateSolution(VectorType **V,
                              const hiflow::la::SeqDenseMatrix< DataType > &H,
                              const std::vector< DataType > &g, int k,
                              VectorType *x) const;

  /// basis of subspace
  std::vector< VectorType * > V_;
  VectorType w_;
  VectorType z_;

  /// max size of the Krylov subspace basis
  int size_basis_;
  int num_init_basis_;
  bool aux_vec_init_;
  bool reuse_basis_;
};

/// @brief GMRES creator class
/// @author Tobias Hahn

template < class LAD > class GMREScreator : public LinearSolverCreator< LAD > {
public:
  LinearSolver< LAD > *params(const PropertyTree &c) {
    GMRES< LAD > *newGMRES = new GMRES< LAD >();
    if (c.contains("Method") && c.contains("SizeBasis")) {
      newGMRES->InitParameter(
          c["SizeBasis"].template get< int >(),
          c["Method"].template get< std::string >().c_str());
    }
    if (c.contains("MaxIterations") && c.contains("AbsTolerance") &&
        c.contains("RelTolerance") && c.contains("DivTolerance")) {
      newGMRES->InitControl(c["MaxIterations"].template get< int >(),
                            c["AbsTolerance"].template get< double >(),
                            c["RelTolerance"].template get< double >(),
                            c["DivTolerance"].template get< double >());
    }
    return newGMRES;
  }
};

template < class LAD, class PreLAD >
void setup_GMRES_solver(GMRES< LAD, PreLAD > &gmres_solver,
                        const PropertyTree &params,
                        NonlinearProblem< LAD > *nonlin) {
  const int max_it = params["MaxIt"].get< int >(1000);
  const int max_size = params["KrylovSize"].get< int >(100);
  const double abs_tol = params["AbsTol"].get< double >(1e-12);
  const double rel_tol = params["RelTol"].get< double >(1e-10);
  const bool use_press_filter = params["UsePressureFilter"].get< bool >(false);

  gmres_solver.InitControl(max_it, abs_tol, rel_tol, 1e6);
  if (params["UsePrecond"].get< bool >(true)) {
    gmres_solver.InitParameter(max_size, "RightPreconditioning");
  } else {
    gmres_solver.InitParameter(max_size, "NoPreconditioning");
  }
  gmres_solver.SetPrintLevel(params["PrintLevel"].get< int >(0));
  gmres_solver.SetReuse(params["Reuse"].get< bool >(true));
  gmres_solver.SetReuseBasis(params["ReuseBasis"].get< bool >(true));

  gmres_solver.SetName(params["Name"].get< std::string >("GMRES"));

  if (use_press_filter && nonlin != nullptr) {
    gmres_solver.SetupNonLinProblem(nonlin);
  }
}

template < class LAD >
void setup_GMRES_solver(GMRES< LAD > &gmres_solver, const PropertyTree &params,
                        NonlinearProblem< LAD > *nonlin) {
  setup_GMRES_solver< LAD, LAD >(gmres_solver, params, nonlin);
}

} // namespace la
} // namespace hiflow

#endif // HIFLOW_LINEARSOLVER_GMRES_H_
