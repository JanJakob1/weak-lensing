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

#include "newton.h"
#include "common/log.h"
#include "common/macros.h"

#include <cassert>
#include <iomanip>
#include <vector>
#include <sstream>
#include <string>

#include "linear_algebra/la_descriptor.h"
#include "linear_algebra/pce_matrix.h"
#include "linear_algebra/block_matrix.h"

#include "damping_strategy.h"
#include "forcing_strategy.h"
#include "nonlinear_problem.h"

#include "assembly/assembly.h"
#include "space/cell_visualization.h"

/// @author Tobias Hahn, Michael Schick

using namespace hiflow::la;

namespace hiflow {

/// Sets up paramaters for initial solution, forcing and damping

template < class LAD >
NonlinearSolverState
Newton< LAD >::InitParameter(NonlinearSolverParameter param) {
  if (param == NewtonDampingStrategyNone || param == NewtonDampingStrategyOwn) {
    this->DampingStrategy_ = param;
  } else if (param == NewtonForcingStrategyConstant ||
             param == NewtonForcingStrategyOwn) {
    this->ForcingStrategy_ = param;
  } else if (param == NewtonInitialSolution0 ||
             param == NewtonInitialSolutionOwn) {
    this->InitialSolution_ = param;
  } else {
    return kNonlinearSolverInitError;
  }

  return kNonlinearSolverSuccess;
}

template < class LAD >
void Newton< LAD >::InitParameter(VectorType *residual, MatrixType *matrix) {
  res_ = residual;
  jac_ = matrix;
  assert(res_ != nullptr);
  assert(jac_ != nullptr);
}

/// Sets up the damping strategy
/// @param dampstrat DampingStrategy

template < class LAD >
void Newton< LAD >::SetDampingStrategy(DampingStrategy< LAD > &dampstrat) {
  this->DampStratObject_ = &dampstrat;
  this->DampingStrategy_ = NewtonDampingStrategyOwn;
}

/// Sets up the forcing strategy
/// @param forcingstrat ForcingStrategy

template < class LAD >
void Newton< LAD >::SetForcingStrategy(ForcingStrategy< LAD > &forcingstrat) {
  this->ForcingStratObject_ = &forcingstrat;
  this->ForcingStrategy_ = NewtonForcingStrategyOwn;
}

/// Get the current forcing term
/// returns zero if no forcing is activated

template < class LAD >
void Newton< LAD >::GetForcingTerm(DataType &forcing) const {
  if (this->ForcingStrategy_ == NewtonForcingStrategyOwn) {
    forcing = this->ForcingStratObject_->GetCurrentForcingTerm();
  } else {
    forcing = 0.;
  }
}

/// Set a forcing term, necessary in combination with damping
/// only if forcing is activated
/// @param forcing DataType

template < class LAD > void Newton< LAD >::SetForcingTerm(DataType forcing) {
  if (this->ForcingStrategy_ == NewtonForcingStrategyOwn) {
    this->ForcingStratObject_->SetForcingTerm(forcing);
  }
}

/// Provides information is forcing is used
/// necessary for damping

template < class LAD > bool Newton< LAD >::Forcing() const {
  return (this->ForcingStrategy_ == NewtonForcingStrategyOwn);
}

/// Updates solution vector using correction and possible damping/
/// forcing strategies that use in turn a right-hand side
/// @param cor correction vector
/// @param rhs right-hand side vector
/// @param sol solution vector

template < class LAD >
NonlinearSolverState
Newton< LAD >::UpdateSolution(const VectorType &cor, const VectorType &rhs,
                              VectorType *sol,
                              VectorSpace< DataType > const *space) {
  if (this->DampingStrategy_ == NewtonDampingStrategyNone) {
    assert(sol->size_local() == cor.size_local());
    assert(sol->size_global() == cor.size_global());
    sol->Axpy(cor, static_cast< DataType >(-1.));
    sol->Update();

    // space provided -> interpolate possible hanging dofs
    if (space != nullptr) 
    {
      if (this->print_level_ >= 3) 
      {
        LOG_INFO("Interpolate hanging dofs before filter ", 1);
      }
      interpolate_constrained_vector< LAD >(*space, *sol);
      sol->Update();
    }

    this->op_->ApplyFilter(*sol);
    sol->Update();
    if (space != nullptr) 
    {
      if (this->print_level_ >= 3) 
      {
        LOG_INFO("Interpolate hanging dofs after filter  ", 1);
      }
      interpolate_constrained_vector< LAD >(*space, *sol);
      sol->Update();
    }

    if (non_const_mode_) 
    {
      this->ComputeResidualNonConst(*sol, rhs, this->res_);
    } 
    else 
    {
      this->ComputeResidual(*sol, rhs, this->res_);
    }

    return kNonlinearSolverSuccess;
  } 
  else if (this->DampingStrategy_ == NewtonDampingStrategyOwn) 
  {
    assert(DampStratObject_ != nullptr);
    DampingState state =
        DampStratObject_->Update(cor, rhs, this->res_, sol, this, space);
    this->residual_ = DampStratObject_->GetResidual();
    this->resids_.push_back(this->residual_);

    if (state != 0) 
    {
      return kNonlinearSolverError;
    }
  } 
  else 
  {
    return kNonlinearSolverError;
  }
  return kNonlinearSolverSuccess;
}

/// Updates solution vector using correction and possible damping/
/// forcing strategies
/// @param cor correction vector
/// @param sol solution vector

template < class LAD >
NonlinearSolverState
Newton< LAD >::UpdateSolution(const VectorType &cor, VectorType *sol,
                              VectorSpace< DataType > const *space) {
  VectorType *rhs = new VectorType();
  rhs->Clear();
  NonlinearSolverState state = this->UpdateSolution(cor, *rhs, sol, space);
  delete rhs;
  return state;
}

/// Returns jacobian matrix J of nonlinear problem F at x
/// @param x point of evaluation
/// @param jacobian at x

template < class LAD >
void Newton< LAD >::ComputeJacobian(const VectorType &x, MatrixType *jacobian) {
  this->op_->EvalGrad(x, jacobian);
}

/// Returns jacobian matrix J of nonlinear problem F at x
/// @param x point of evaluation
/// @param jacobian at x

template < class LAD >
void Newton< LAD >::ComputeJacobianNonConst(VectorType &x,
                                            MatrixType *jacobian) {
  this->op_->EvalGradNonConst(x, jacobian);
}

/// Solves linear problem J*c=r
/// If Forcing is activated, then the system
/// is solved in an inexact way with
/// relative tolerance determined by forcing terms
/// @param jacobian jacobian matrix J
/// @param residual residual vector r
/// @param correction correction vector c

template < class LAD >
LinearSolverState Newton< LAD >::SolveJacobian(const MatrixType &jacobian,
                                               const VectorType &residual,
                                               VectorType *correction) {
  assert(correction != nullptr);

  // Reset correction if needed
  if ((residual.size_local() != correction->size_local()) ||
      (residual.size_global() != correction->size_global())) {
    correction->CloneFromWithoutContent(residual);
  }

  // start vector
  correction->Zeros();
  //      correction->Update( );

  if (this->ForcingStrategy_ == NewtonForcingStrategyOwn) {
    assert(this->ForcingStratObject_ != nullptr);
    if (this->print_level_ >= 3) {
      LOG_INFO("Forcing term                      ",
               this->ForcingStratObject_->GetCurrentForcingTerm());
    }
    this->linsolve_->SetRelativeTolerance(
        this->ForcingStratObject_->GetCurrentForcingTerm());
  }

  // pass on information object
  if (this->info_ != nullptr) {
    std::stringstream lin_str;
    lin_str << "linsolve" << std::setw(1 + std::log10(this->control().maxits()))
            << std::setfill('0') << this->iter_;
    this->info_->add(lin_str.str());
    this->linsolve_->SetInfo(this->info_->get_child(lin_str.str()));
  }

  // solve
  LinearSolverState state = this->linsolve_->Solve(residual, correction);

  correction->Update();

  return state; // solve jacobian
}

/// Computes residual vector F(sol)-rhs for non-linear problem F with
/// right-hand side rhs
/// @param sol solution vector
/// @param rhs right-hand side vector
/// @param res residual vector

template < class LAD >
void Newton< LAD >::ComputeResidual(const VectorType &sol,
                                    const VectorType &rhs, VectorType *res) 
{
  assert(res != nullptr);
  // Reset residual if needed
  if ((res->size_local() != sol.size_local()) ||
      (res->size_global() != sol.size_global())) 
  {
    res->CloneFromWithoutContent(sol);
    res->Zeros();
  }

  // Compute residual
  this->op_->EvalFunc(sol, res);
  if ((rhs.size_local() == res->size_local()) &&
      (rhs.size_global() == res->size_global())) 
  {
    res->Axpy(rhs, static_cast< DataType >(-1.));
  }

  // res->Update( );

  // Compute new residual norm
  this->residual_ = res->Norm2();

  // Store it in vector
  this->resids_.push_back(this->residual_);
}

/// Computes residual vector F(sol)-rhs for non-linear problem F with
/// right-hand side rhs
/// @param sol solution vector
/// @param rhs right-hand side vector
/// @param res residual vector

template < class LAD >
void Newton< LAD >::ComputeResidualNonConst(VectorType &sol, const VectorType &rhs, VectorType *res) 
{
  assert(res != nullptr);
  // Reset residual if needed
  if ((res->size_local() != sol.size_local()) ||
      (res->size_global() != sol.size_global())) 
  {
    res->CloneFromWithoutContent(sol);
    res->Zeros();
  }

  // Compute residual
  this->op_->EvalFuncNonConst(sol, res);
  if ((rhs.size_local() == res->size_local()) &&
      (rhs.size_global() == res->size_global())) {
    res->Axpy(rhs, static_cast< DataType >(-1.));
  }

  // res->Update( );

  // Compute new residual norm
  this->residual_ = res->Norm2();

  // Store it in vector
  this->resids_.push_back(this->residual_);
}

/// Computes residual vector F(sol) for non-linear problem F
/// @param sol solution vector
/// @param res residual vector

template < class LAD >
void Newton< LAD >::ComputeResidual(const VectorType &sol, VectorType *res) 
{
  VectorType *rhs = new VectorType();
  rhs->Clear();
  this->ComputeResidual(sol, *rhs, res);
  delete rhs;
}

/// Solves F(x)=y
/// @param y right hand side vectorNewtonDampingStrategyArmijo
/// @param x solution vector
/// @return status if solver succeeded

template < class LAD >
NonlinearSolverState
Newton< LAD >::Solve(const VectorType &rhs, VectorType *x, VectorSpace< DataType > const *space) 
{
  assert(this->res_ != nullptr);
  assert(this->jac_ != nullptr);
  assert(this->op_ != nullptr);
  assert(this->linsolve_ != nullptr);
  assert(x != nullptr);

  Timer timer;
  if (this->info_ != nullptr) {
    timer.reset();
    timer.start();
  }

  // Init
  LinearSolverState LinSolState = kSolverSuccess;
  IterateControl::State conv = IterateControl::kIterate;

  this->res_->Clear();
  this->res_->CloneFromWithoutContent(rhs);

  VectorType *cor = new VectorType();
  cor->Clear();
  cor->CloneFromWithoutContent(rhs);

  VectorType *sol = new VectorType();
  sol->Clear();
  if (InitialSolution_ == NewtonInitialSolutionOwn) 
  {
    sol->CloneFrom(*x);
  } 
  else if (InitialSolution_ == NewtonInitialSolution0) 
  {
    sol->CloneFromWithoutContent(rhs);
    sol->Zeros();
  } 
  else 
  {
    return kNonlinearSolverInitError;
  }

  // Step 0
  this->iter_ = 0;
  this->op_->Reinit();

  sol->Update();

  if (space != nullptr) 
  {
    if (this->print_level_ >= 3) 
    {
      LOG_INFO("Interpolate hanging dofs", 1);
    }
    interpolate_constrained_vector< LAD >(*space, *sol);
    sol->Update();
  }

  if (non_const_mode_) 
  {
    this->ComputeResidualNonConst(*sol, rhs, this->res_);
  } 
  else 
  {
    this->ComputeResidual(*sol, rhs, this->res_);
  }

  conv = this->control().Check(this->iter(), this->GetResidual());
  if (this->ForcingStrategy_ == NewtonForcingStrategyOwn) 
  {
    assert(this->ForcingStratObject_ != nullptr);
    this->ForcingStratObject_->SetResidualNorm(this->res_->Norm2());
  }

  DataType lin_solver_rhs_norm = this->GetResidual();
  if (this->print_level_ >= 1) 
  {
    LOG_INFO("Initial res norm", this->GetResidual());
  }
  if (this->ForcingStrategy_ == NewtonForcingStrategyOwn) 
  {
    if (this->print_level_ >= 1) 
    {
      LOG_INFO("Initial forcing term", this->ForcingStratObject_->GetCurrentForcingTerm());
    }
  }

  while (conv == IterateControl::kIterate) 
  {
    // NextStep
    this->iter_++;
    Timer timer;
    timer.start();

    if (non_const_mode_) 
    {
      this->ComputeJacobianNonConst(*sol, this->jac_);
    } 
    else 
    {
      this->ComputeJacobian(*sol, this->jac_);
    }
    timer.stop();
    DataType asm_time = timer.get_duration();
    timer.reset();

    lin_solver_rhs_norm = this->GetResidual();

    timer.start();
    LinSolState = this->SolveJacobian(*this->jac_, *this->res_, cor);

    timer.stop();
    DataType solve_time = timer.get_duration();
    timer.reset();
    timer.start();

    if (LinSolState == kSolverError) 
    {
      break;
    }
    this->op_->ApplyCorrectionFilter(*cor);

    this->UpdateSolution(*cor, rhs, sol, space);
    timer.stop();
    DataType update_time = timer.get_duration();

    conv = this->control().Check(this->iter(), this->GetResidual());

    int lin_iter = this->linsolve_->iter();
    DataType lin_res = this->linsolve_->res();

    if (this->print_level_ >= 1) 
    {
      LOG_INFO("", "==================================");
      LOG_INFO("Iteration", this->iter_);
      LOG_INFO("Abs res norm", this->GetResidual());
    }
    if (this->print_level_ >= 2) 
    {
      LOG_INFO("LinSolver iter", lin_iter);
      LOG_INFO("LinSolver abs res", lin_res);
      LOG_INFO("LinSolver rel res", lin_res / lin_solver_rhs_norm);
    }

    if (this->ForcingStrategy_ == NewtonForcingStrategyOwn) 
    {
      assert(this->ForcingStratObject_ != nullptr);
      this->ForcingStratObject_->ComputeForcingTerm(this->res_->Norm2(), this->linsolve_->res());
      if (this->print_level_ >= 2) 
      {
        LOG_INFO("New forcing term", this->ForcingStratObject_->GetCurrentForcingTerm());
      }
    }

    if (sol->my_rank() == 0) 
    {
      this->WriteStatistics(asm_time, solve_time, update_time, lin_solver_rhs_norm);
    }
    this->asm_time_ += asm_time;
    this->solve_time_ += solve_time;
    this->update_time_ += update_time;
  }
  delete cor;

  if (this->print_level_ >= 1) 
  {
    LOG_INFO("Final res norm", this->GetResidual());
  }
  if (this->info_ != nullptr) 
  {
    timer.stop();
    this->info_->add("iter", this->iter());
    this->info_->add("time", timer.get_duration());
  }

  if (LinSolState == kSolverError) 
  {
    if (this->force_return_sol_) 
    {
      x->CopyFrom(*sol);
      x->Update();
    }
    delete sol;
    return kNonlinearSolverError;
  } 
  else if (conv == IterateControl::kFailureDivergenceTol ||
             conv == IterateControl::kFailureMaxitsExceeded) 
  {
    if (this->force_return_sol_) 
    {
      x->CopyFrom(*sol);
      x->Update();
    }
    delete sol;
    return kNonlinearSolverExceeded;
  } 
  else 
  {
    x->CopyFrom(*sol);
    x->Update();
    delete sol;
    return kNonlinearSolverSuccess;
  }
}

/// Solves F(x)=0
/// @param x solution vector
/// @return status if solver succeeded

template < class LAD >
NonlinearSolverState
Newton< LAD >::Solve(VectorType *x, VectorSpace< DataType > const *space) {
  VectorType *rhs = new VectorType();
  rhs->Clear();
  rhs->CloneFromWithoutContent(*x);
  rhs->Zeros();
  NonlinearSolverState state = this->Solve(*rhs, x, space);
  delete rhs;
  return state;
}

template < class LAD >
void Newton< LAD >::WriteStatistics(DataType asm_time, DataType solve_time,
                                    DataType update_time,
                                    DataType lin_solver_rhs_norm) {
  if (this->filename_statistics_.empty()) {
    return;
  }
  DataType lin_iter = this->linsolve_->iter();
  DataType lin_res = this->linsolve_->res();
  DataType newton_iter = this->iter_;
  int resids_size = this->resids_.size();
  DataType newton_res = this->resids_[resids_size - 1];

  std::string path = this->filename_statistics_;
  std::ofstream out;
  out.open(path.c_str(), std::ios::out | std::ios::app);
  out.precision(6);
  out << std::scientific;
  out << newton_iter << " " << newton_res << " " << lin_iter << " " << lin_res
      << " " << lin_res / lin_solver_rhs_norm << " " << asm_time << " "
      << solve_time << " " << update_time << " "
      << "\n";
  out.close();
}

template < class LAD >
void Newton< LAD >::visualize_solution (VectorType &sol, 
                                        const VectorSpace< DataType > *space,
                                        std::string const &prefix, 
                                        int iter) const
{ 
  if (space == nullptr)
  {
    return;
  }
  
  LOG_INFO("Visualize Newton iteration", iter);

  std::stringstream input;

  if (iter < 10) {
    input << prefix << ".000" << iter;
  } 
  else if (iter < 100) {
    input << prefix << ".00" << iter;
  } else if (iter < 1000) {
    input << prefix << ".0" << iter;
  } else {
    input << prefix << "." << iter;
  }

  int num_partitions = space->dof().get_nb_subdomains();
  int num_variables = space->get_nb_var();
  if (num_partitions > 1) 
  {
    input << ".pvtu";
  } 
  else {
    input << ".vtu";
  }

  std::string visu_filename = input.str();

  ParallelCellVisualization< DataType > visu(*space, 1, space->dof().get_mpi_comm(), 0);
  
  sol.Update();
  for (int v=0; v<num_variables; ++v)
  {
    std::stringstream var_name;
    var_name << "var_" << v;
    visu.visualize(FeEvalOnCellScalar< LAD >(*space, sol, v), var_name.str());
  }
  
  // parallel statistics
  std::vector< DataType > remote_index(space->meshPtr()->num_entities(space->meshPtr()->tdim()), 0);
  std::vector< DataType > sub_domain(space->meshPtr()->num_entities(space->meshPtr()->tdim()), 0);

  AttributePtr sub = space->meshPtr()->get_attribute("_sub_domain_", space->meshPtr()->tdim());
  AttributePtr remote = space->meshPtr()->get_attribute("_remote_index_", space->meshPtr()->tdim());
  for (mesh::EntityIterator it = space->meshPtr()->begin(space->meshPtr()->tdim());
       it != space->meshPtr()->end(space->meshPtr()->tdim()); ++it) 
  {
    remote_index.at(it->index()) = remote->get_int_value(it->index());
    sub_domain.at(it->index()) = sub->get_int_value(it->index());
  }
  visu.visualize_cell_data(remote_index, "_remote_index_");
  visu.visualize_cell_data(sub_domain, "_sub_domain_");

  // write file
  visu.write(visu_filename);
}  
 
/// Plain constructor that does not initialize the solver.

template < class LAD > Newton< LAD >::Newton() {
  this->InitialSolution_ = NewtonInitialSolution0;
  this->DampingStrategy_ = NewtonDampingStrategyNone;
  this->ForcingStrategy_ = NewtonForcingStrategyConstant;
  this->non_const_mode_ = false;

  this->asm_time_ = 0.;
  this->solve_time_ = 0.;
  this->update_time_ = 0.;
  this->filename_statistics_ = "";
}

/// Standard constructor requiring pointers to user reserved space
/// for residual vector and jacobian. Sets no forcing and no damping and
/// uses initial solution zero.

template < class LAD >
Newton< LAD >::Newton(VectorType *residual, MatrixType *matrix)
    : res_(residual), jac_(matrix) {
  assert(res_ != nullptr);
  assert(jac_ != nullptr);

  this->InitialSolution_ = NewtonInitialSolution0;
  this->DampingStrategy_ = NewtonDampingStrategyNone;
  this->ForcingStrategy_ = NewtonForcingStrategyConstant;
  this->non_const_mode_ = false;

  this->asm_time_ = 0.;
  this->solve_time_ = 0.;
  this->update_time_ = 0.;
  this->filename_statistics_ = "";
  this->force_return_sol_ = false;
}

/// Standard destructor

template < class LAD > Newton< LAD >::~Newton() {
  this->res_ = nullptr;
  this->jac_ = nullptr;
  this->linsolve_ = nullptr;
  this->op_ = nullptr;
  this->DampStratObject_ = nullptr;
  this->ForcingStratObject_ = nullptr;
  this->info_ = nullptr;
  this->filename_statistics_ = "";
  this->force_return_sol_ = false;
}

/// template instantiation
template class Newton< LADescriptorCoupledD >;
template class Newton< LADescriptorCoupledS >;
template class Newton< LADescriptorBlock< LADescriptorCoupledD > >;
template class Newton< LADescriptorBlock< LADescriptorCoupledS > >;
template class Newton< LADescriptorPCE< LADescriptorCoupledD > >;
template class Newton< LADescriptorPCE< LADescriptorCoupledS > >;
template class Newton< LADescriptorPolynomialChaosD >;

#ifdef WITH_HYPRE
template class Newton< LADescriptorHypreD >;
template class Newton< LADescriptorPCE< LADescriptorHypreD > >;
template class Newton< LADescriptorBlock< LADescriptorHypreD > >;
#endif

} // namespace hiflow
