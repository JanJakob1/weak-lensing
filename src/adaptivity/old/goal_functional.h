// Copyright (C) 2011-2017 Vincent Heuveline
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

#ifndef HIFLOW_ADAPTIVITY_GOAL_FUNCTIONAL_H
#define HIFLOW_ADAPTIVITY_GOAL_FUNCTIONAL_H

#include "assembly/function_values.h"
#include "common/vector_algebra.h"
#include "linear_algebra/la_descriptor.h"
#include <iostream>
#include <string>
#include <cmath>
#include <iostream>
#include <limits>

///
/// \file goal_functional.h Abstract structure for Goal-functionals.
/// \author Martin Baumann, Philipp Gerstner
///
/// For a goal-functional \f$J\f$ of the following form
///   \f[
///        J(u):=\int_0^T (j_1(t,u(t)),\varphi)_\Omega dt +
///                       (j_2(T,u(T)),\varphi)_\Omega.
///   \f]
/// The part \f$j_1\f$ is denoted by the function j_force_type and its related
/// parameters are stored in the struct force_type_parameters.
///

using hiflow::Vec;
using namespace hiflow::la;

namespace hiflow {

template < int DIM, class DataType > struct ParametersForceType {
  Vec< DIM, DataType > x;
  DataType absolute_time;
  int var;
  DataType relative_time; // -1 -> t=t_{i-1}, 0 -> t=t_i, 1->t=t_{i+1}
  std::vector< DataType > solP_prev;                  // primal solution t_{i-1}
  std::vector< DataType > solP;                       // primal solution t_i
  std::vector< DataType > solP_next;                  // primal solution t_i+1
  std::vector< Vec< DIM, DataType > > grad_solP_prev; // primal solution t_{i-1}
  std::vector< Vec< DIM, DataType > > grad_solP;      // primal solution t_i
  std::vector< Vec< DIM, DataType > > grad_solP_next; // primal solution t_i+1
  DataType phi;                  // test function for variable 'var'
  Vec< DIM, DataType > grad_phi; // gradient of test function for variable 'var'
};

template < int DIM, class DataType > struct ParametersFinalType {
  Vec< DIM, DataType > x;
  int var;
  DataType absolute_time;
  std::vector< DataType > solP;                  // primal solution
  std::vector< Vec< DIM, DataType > > grad_solP; // primal solution
  DataType phi;                  // test function for variable 'var'
  Vec< DIM, DataType > grad_phi; // gradient of test function for variable 'var'
};

template < int DIM, class DataType > struct ParametersEvalType {
  Vec< DIM, DataType > x;
  DataType absolute_time;
  int var;
  std::vector< DataType > solP;                    // primal solution t_i
  std::vector< Vec< DIM, DataType > > grad_solP;   // primal solution t_i
  std::vector< DataType > d_solP;                  // primal solution t_i
  std::vector< Vec< DIM, DataType > > d_grad_solP; // primal solution t_i
};

template < int DIM, class DataType > 
class GoalFunctional {
public:
  GoalFunctional();

  virtual ~GoalFunctional() { ; }

  /// \brief Force-type contribution of a goal-functional
  /// Only the part \f$(j_1(t,u(t)),varphi)\f$ should be implemented, excluding
  /// parts related to quadrature in space or in time.
  virtual DataType
  j_force_type(ParametersForceType< DIM, DataType > p) const = 0;

  /// \brief final-type contribution of a goal-functional (related to \f$t=T\f$)
  /// Only the part \f$(j_1(t,u(t)),varphi)\f$ should be implemented, excluding
  /// parts related to quadrature in space or in time.
  virtual DataType
  j_final_type(ParametersFinalType< DIM, DataType > p) const = 0;

  virtual DataType
  j_force_eval(ParametersEvalType< DIM, DataType > p) const = 0;

  virtual DataType
  j_final_eval(ParametersEvalType< DIM, DataType > p) const = 0;

  virtual DataType
  j_force_eval_deriv(ParametersEvalType< DIM, DataType > p) const {
    return 0.;
  };

  virtual DataType
  j_final_eval_deriv(ParametersEvalType< DIM, DataType > p) const {
    return 0.;
  };

  bool &force_type_active() { return force_type_active_; }

  bool &final_type_active() { return final_type_active_; }

  void set_active_parts(bool active_force, bool active_ic) {
    this->force_type_active_ = active_force;
    this->final_type_active_ = active_ic;
  }

  void set_cyl_box(DataType t_min, DataType t_max, DataType z_min,
                   DataType z_max, DataType r_min, DataType r_max,
                   DataType phi_min, DataType phi_max) {
    this->t_min_ = t_min;
    this->t_max_ = t_max;
    this->z_min_ = z_min;
    this->z_max_ = z_max;
    this->r_min_ = r_min;
    this->r_max_ = r_max;
    this->phi_min_ = phi_min;
    this->phi_max_ = phi_max;
    this->co_system_ = 1;
  }

  void set_cart_box(DataType t_min, DataType t_max, DataType x_min,
                    DataType x_max, DataType y_min, DataType y_max,
                    DataType z_min, DataType z_max) {
    this->t_min_ = t_min;
    this->t_max_ = t_max;
    this->z_min_ = z_min;
    this->z_max_ = z_max;
    this->y_min_ = y_min;
    this->y_max_ = y_max;
    this->x_min_ = x_min;
    this->x_max_ = x_max;
    this->co_system_ = 0;
  }

  void set_scale_factor(DataType scale) { this->scale_ = scale; }

  int GetIntegralType() { return this->integral_type_; }

  void use_cart_co_system() { this->co_system_ = 0; }

  void use_cyl_co_system() { this->co_system_ = 1; }

protected:
  bool check_force_type(const ParametersForceType< DIM, DataType > &p) const;
  bool check_force_eval(const ParametersEvalType< DIM, DataType > &p) const;
  bool check_final_type(const ParametersFinalType< DIM, DataType > &p) const;
  bool check_final_eval(const ParametersEvalType< DIM, DataType > &p) const;

  bool is_in_box(Vec< DIM, DataType > x) const;

  void extract_solP_force(const ParametersForceType< DIM, DataType > &p,
                          std::vector< DataType > &solP) const;

  void
  extract_grad_solP_force(const ParametersForceType< DIM, DataType > &p,
                          std::vector< Vec< DIM, DataType > > &grad_solP) const;

  std::string name_;
  bool force_type_active_;
  bool final_type_active_;

  DataType x_min_;
  DataType x_max_;
  DataType y_min_;
  DataType y_max_;
  DataType z_min_;
  DataType z_max_;
  DataType r_min_;
  DataType r_max_;
  DataType phi_min_;
  DataType phi_max_;
  DataType t_min_;
  DataType t_max_;
  DataType scale_;

  int co_system_;
  int integral_type_;
};

// **************************************************************************
// Specific veriable in specifc subdomain during specific time interval
// **************************************************************************

template < int DIM, class DataType >
class GoalFunctionalVariableSubDomain : public GoalFunctional< DIM, DataType > {
public:
  GoalFunctionalVariableSubDomain();

  void set_variable(int var) { this->var_ = var; }

  DataType j_force_type(ParametersForceType< DIM, DataType > p) const;
  DataType j_final_type(ParametersFinalType< DIM, DataType > p) const;

  DataType j_force_eval(ParametersEvalType< DIM, DataType > p) const;
  DataType j_final_eval(ParametersEvalType< DIM, DataType > p) const;
  DataType j_force_eval_deriv(ParametersEvalType< DIM, DataType > p) const;
  DataType j_final_eval_deriv(ParametersEvalType< DIM, DataType > p) const;

protected:
  int var_;
};

// **************************************************************************
// Specific derivative of specific veriable in specifc subdomain during specific
// time interval
// **************************************************************************

template < int DIM, class DataType >
class GoalFunctionalDerivativeSubDomain
    : public GoalFunctional< DIM, DataType > {
public:
  GoalFunctionalDerivativeSubDomain();

  void set_variable(int var) { this->var_ = var; }

  void set_derivative(int deriv) {
    assert(deriv < DIM);
    this->deriv_ = deriv;
  }

  void set_squared(bool flag) { this->squared_ = flag; }

  DataType j_force_type(ParametersForceType< DIM, DataType > p) const;
  DataType j_final_type(ParametersFinalType< DIM, DataType > p) const;

  DataType j_force_eval(ParametersEvalType< DIM, DataType > p) const;
  DataType j_final_eval(ParametersEvalType< DIM, DataType > p) const;
  DataType j_force_eval_deriv(ParametersEvalType< DIM, DataType > p) const;
  DataType j_final_eval_deriv(ParametersEvalType< DIM, DataType > p) const;

protected:
  bool squared_;
  int var_;
  int deriv_;
};


////////////////////////////////////////////////////////////////////////////////
//////////////////////////// Implementation ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template < int DIM, class DataType >
GoalFunctional< DIM, DataType >::GoalFunctional() {
  force_type_active_ = false;
  final_type_active_ = false;
  name_ = "NOT SET";
  this->scale_ = 1.;
  this->co_system_ = 0;
  this->t_min_ = 0.;
  this->t_max_ = 99999.;
}

template < int DIM, class DataType >
bool GoalFunctional< DIM, DataType >::is_in_box(Vec< DIM, DataType > x) const {
  switch (this->co_system_) {
  case 0:
    if (x[0] < this->x_min_ || x[0] > this->x_max_)
      return false;

    if (DIM > 1) {
      if (x[1] < this->y_min_ || x[1] > this->y_max_)
        return false;

      if (DIM > 2) {
        if (x[2] < this->z_min_ || x[2] > this->z_max_)
          return false;
      }
    }
    break;
  case 1:
    if (x[0] < this->phi_min_ || x[0] > this->phi_max_)
      return false;

    if (x[1] < this->r_min_ || x[1] > this->r_max_)
      return false;

    if (DIM > 2) {
      if (x[2] < this->z_min_ || x[2] > this->z_max_)
        return 0.;
    }
    break;
  }
  return true;
}

template < int DIM, class DataType >
void GoalFunctional< DIM, DataType >::extract_solP_force(
    const ParametersForceType< DIM, DataType > &p,
    std::vector< DataType > &solP) const {
  if (p.relative_time == -1) {
    solP = p.solP_prev;
  } else if (p.relative_time == -0.5) {
    solP = p.solP_prev;
    for (int l = 0; l < DIM; ++l) {
      solP[l] += p.solP[l];
      solP[l] *= 0.5;
    }
  } else if (p.relative_time == 0) {
    solP = p.solP;
  } else if (p.relative_time == 0.5) {
    solP = p.solP_next;
    for (int l = 0; l < DIM; ++l) {
      solP[l] += p.solP[l];
      solP[l] *= 0.5;
    }
  }

  else if (p.relative_time == 1) {
    solP = p.solP_next;
  }
}

template < int DIM, class DataType >
void GoalFunctional< DIM, DataType >::extract_grad_solP_force(
    const ParametersForceType< DIM, DataType > &p,
    std::vector< Vec< DIM, DataType > > &grad_solP) const {
  if (p.relative_time == -1) {
    grad_solP = p.grad_solP_prev;
  } else if (p.relative_time == -0.5) {
    grad_solP = p.grad_solP_prev;
    for (int l = 0; l < DIM; ++l) {
      for (int d = 0; d < DIM; ++d) {
        grad_solP[l][d] += p.grad_solP[l][d];
        grad_solP[l][d] *= 0.5;
      }
    }
  } else if (p.relative_time == 0) {
    grad_solP = p.grad_solP;
  } else if (p.relative_time == 0.5) {
    grad_solP = p.grad_solP_next;
    for (int l = 0; l < DIM; ++l) {
      for (int d = 0; d < DIM; ++d) {
        grad_solP[l][d] += p.grad_solP[l][d];
        grad_solP[l][d] *= 0.5;
      }
    }
  } else if (p.relative_time == 1) {
    grad_solP = p.grad_solP_next;
  }
}

template < int DIM, class DataType >
bool GoalFunctional< DIM, DataType >::check_force_type(
    const ParametersForceType< DIM, DataType > &p) const {
  if (!this->force_type_active_) {
    return false;
  }
  if (!this->is_in_box(p.x)) {
    return false;
  }
  if (p.absolute_time < this->t_min_ || p.absolute_time > this->t_max_) {
    return false;
  }
  return true;
}

template < int DIM, class DataType >
bool GoalFunctional< DIM, DataType >::check_force_eval(
    const ParametersEvalType< DIM, DataType > &p) const {
  if (!this->force_type_active_) {
    return false;
  }
  if (!this->is_in_box(p.x)) {
    return false;
  }
  if (p.absolute_time < this->t_min_ || p.absolute_time > this->t_max_) {
    return false;
  }
  return true;
}

template < int DIM, class DataType >
bool GoalFunctional< DIM, DataType >::check_final_type(
    const ParametersFinalType< DIM, DataType > &p) const {
  if (!this->final_type_active_) {
    return false;
  }
  if (!this->is_in_box(p.x)) {
    return false;
  }
  return true;
}

template < int DIM, class DataType >
bool GoalFunctional< DIM, DataType >::check_final_eval(
    const ParametersEvalType< DIM, DataType > &p) const {
  if (!this->final_type_active_) {
    return false;
  }
  if (!this->is_in_box(p.x)) {
    return false;
  }
  return true;
}

template class GoalFunctional< 2, double >;
template class GoalFunctional< 3, double >;

// **************************************************************************
// Certain Variable in specifc subdomain of cylindrical geometry inside specific
// time interval
// **************************************************************************

template < int DIM, class DataType >
GoalFunctionalVariableSubDomain< DIM,
                                 DataType >::GoalFunctionalVariableSubDomain() {
  this->name_ = "VariableSubDomain";

  this->force_type_active_ = true;
  this->final_type_active_ = true;
  this->integral_type_ = 1;
  this->var_ = 0;
}

template < int DIM, class DataType >
DataType GoalFunctionalVariableSubDomain< DIM, DataType >::j_force_type(
    ParametersForceType< DIM, DataType > p) const {
  if (!this->check_force_type(p)) {
    return 0.;
  }
  if (p.var != var_) {
    return 0.;
  }
  return p.phi * this->scale_;
}

template < int DIM, class DataType >
DataType GoalFunctionalVariableSubDomain< DIM, DataType >::j_force_eval(
    ParametersEvalType< DIM, DataType > p) const {
  if (!this->check_force_eval(p)) {
    return 0.;
  }
  return p.solP[var_] * this->scale_;
}

template < int DIM, class DataType >
DataType GoalFunctionalVariableSubDomain< DIM, DataType >::j_final_type(
    ParametersFinalType< DIM, DataType > p) const {
  if (!this->check_final_type(p)) {
    return 0.;
  }
  if (p.var != var_) {
    return 0.;
  }
  return p.phi * this->scale_;
}

template < int DIM, class DataType >
DataType GoalFunctionalVariableSubDomain< DIM, DataType >::j_final_eval(
    ParametersEvalType< DIM, DataType > p) const {
  if (!this->final_type_active_) {
    return 0.;
  }
  if (!this->is_in_box(p.x)) {
    return 0.;
  }
  return p.solP[var_] * this->scale_;
}

template < int DIM, class DataType >
DataType GoalFunctionalVariableSubDomain< DIM, DataType >::j_force_eval_deriv(
    ParametersEvalType< DIM, DataType > p) const {
  if (!this->check_force_eval(p)) {
    return 0.;
  }
  if (p.var != var_) {
    return 0.;
  }
  return p.d_solP[var_] * this->scale_;
}

template < int DIM, class DataType >
DataType GoalFunctionalVariableSubDomain< DIM, DataType >::j_final_eval_deriv(
    ParametersEvalType< DIM, DataType > p) const {
  if (!this->check_final_eval(p)) {
    return 0.;
  }
  if (p.var != var_) {
    return 0.;
  }
  return p.d_solP[var_] * this->scale_;
}

template class GoalFunctionalVariableSubDomain< 2, double >;
template class GoalFunctionalVariableSubDomain< 3, double >;

// **************************************************************************
// Certain Variable in specifc subdomain of cylindrical geometry inside specific
// time interval
// **************************************************************************

template < int DIM, class DataType >
GoalFunctionalDerivativeSubDomain<
    DIM, DataType >::GoalFunctionalDerivativeSubDomain() {
  this->name_ = "VariableSubDomain";

  this->force_type_active_ = true;
  this->final_type_active_ = true;
  this->integral_type_ = 1;
  this->var_ = 0;
  this->deriv_ = 0;
  this->squared_ = false;
}

template < int DIM, class DataType >
DataType GoalFunctionalDerivativeSubDomain< DIM, DataType >::j_force_type(
    ParametersForceType< DIM, DataType > p) const {
  if (!this->check_force_type(p)) {
    return 0.;
  }
  if (p.var != var_) {
    return 0.;
  }

  std::vector< Vec< DIM, DataType > > grad_solP;
  this->extract_grad_solP_force(p, grad_solP);

  if (this->co_system_ == 0) {
    if (!this->squared_) {
      return p.grad_phi[deriv_] * this->scale_;
    } else {
      return this->scale_ * 2.0 * grad_solP[this->var_][this->deriv_] *
             p.grad_phi[deriv_];
    }
  } else {
    const DataType r = p.x[1];
    if (!this->squared_) {
      const DataType inv_r_comp[3] = {1. / r, 1., 1.};
      return inv_r_comp[deriv_] * p.grad_phi[deriv_] * this->scale_;
    } else {
      const DataType rr = r * r;
      const DataType inv_rr_comp[3] = {1. / rr, 1., 1.};
      return this->scale_ * 2.0 * inv_rr_comp[this->deriv_] *
             grad_solP[this->var_][this->deriv_] * p.grad_phi[deriv_];
    }
  }
}

template < int DIM, class DataType >
DataType GoalFunctionalDerivativeSubDomain< DIM, DataType >::j_force_eval(
    ParametersEvalType< DIM, DataType > p) const {
  if (!this->check_force_eval(p)) {
    return 0.;
  }
  if (this->co_system_ == 0) {
    if (!this->squared_) {
      return p.grad_solP[var_][deriv_] * this->scale_;
    } else {
      return this->scale_ * p.grad_solP[var_][deriv_] *
             p.grad_solP[var_][deriv_];
    }
  } else {
    const DataType r = p.x[1];
    if (!this->squared_) {
      const DataType inv_r_comp[3] = {1. / r, 1., 1.};
      return inv_r_comp[deriv_] * p.grad_solP[var_][deriv_] * this->scale_;
    } else {
      const DataType rr = r * r;
      const DataType inv_rr_comp[3] = {1. / rr, 1., 1.};
      return inv_rr_comp[deriv_] * p.grad_solP[var_][deriv_] *
             p.grad_solP[var_][deriv_] * this->scale_;
    }
  }
}

template < int DIM, class DataType >
DataType GoalFunctionalDerivativeSubDomain< DIM, DataType >::j_final_type(
    ParametersFinalType< DIM, DataType > p) const {
  if (!this->check_final_type(p)) {
    return 0.;
  }
  if (p.var != var_) {
    return 0.;
  }
  if (this->co_system_ == 0) {
    if (!this->squared_) {
      return p.grad_phi[deriv_] * this->scale_;
    } else {
      return this->scale_ * 2.0 * p.grad_solP[this->var_][this->deriv_] *
             p.grad_phi[deriv_];
    }
  } else {
    const DataType r = p.x[1];
    if (!this->squared_) {
      const DataType inv_r_comp[3] = {1. / r, 1., 1.};
      return inv_r_comp[deriv_] * p.grad_phi[deriv_] * this->scale_;
    } else {
      const DataType rr = r * r;
      const DataType inv_rr_comp[3] = {1. / rr, 1., 1.};
      return this->scale_ * 2.0 * inv_rr_comp[deriv_] *
             p.grad_solP[this->var_][this->deriv_] * p.grad_phi[deriv_];
    }
  }
}

template < int DIM, class DataType >
DataType GoalFunctionalDerivativeSubDomain< DIM, DataType >::j_final_eval(
    ParametersEvalType< DIM, DataType > p) const {
  if (!this->check_final_eval(p)) {
    return 0.;
  }
  if (this->co_system_ == 0) {
    if (!this->squared_) {
      return p.grad_solP[var_][deriv_] * this->scale_;
    } else {
      return p.grad_solP[var_][deriv_] * p.grad_solP[var_][deriv_] *
             this->scale_;
    }
  } else {
    const DataType r = p.x[1];
    if (!this->squared_) {
      const DataType inv_r_comp[3] = {1. / r, 1., 1.};
      return inv_r_comp[deriv_] * p.grad_solP[var_][deriv_] * this->scale_;
    } else {
      const DataType rr = r * r;
      const DataType inv_rr_comp[3] = {1. / rr, 1., 1.};
      return inv_rr_comp[deriv_] * p.grad_solP[var_][deriv_] *
             p.grad_solP[var_][deriv_] * this->scale_;
    }
  }
}

template < int DIM, class DataType >
DataType GoalFunctionalDerivativeSubDomain< DIM, DataType >::j_force_eval_deriv(
    ParametersEvalType< DIM, DataType > p) const {
  if (!this->check_force_eval(p)) {
    return 0.;
  }

  if (p.var != var_) {
    return 0.;
  }

  std::vector< Vec< DIM, DataType > > grad_solP = p.grad_solP;

  if (this->co_system_ == 0) {
    if (!this->squared_) {
      return p.d_grad_solP[var_][deriv_] * this->scale_;
    } else {
      return this->scale_ * 2.0 * grad_solP[var_][deriv_] *
             p.d_grad_solP[var_][deriv_];
    }
  } else {
    const DataType r = p.x[1];
    if (!this->squared_) {
      const DataType inv_r_comp[3] = {1. / r, 1., 1.};
      return inv_r_comp[deriv_] * p.d_grad_solP[var_][deriv_] * this->scale_;
    } else {
      const DataType rr = r * r;
      const DataType inv_rr_comp[3] = {1. / rr, 1., 1.};
      return this->scale_ * 2.0 * inv_rr_comp[2] * grad_solP[var_][deriv_] *
             p.d_grad_solP[var_][deriv_];
    }
  }
}

template < int DIM, class DataType >
DataType GoalFunctionalDerivativeSubDomain< DIM, DataType >::j_final_eval_deriv(
    ParametersEvalType< DIM, DataType > p) const {
  if (!this->check_final_eval(p)) {
    return 0.;
  }

  if (p.var != var_) {
    return 0.;
  }
  if (this->co_system_ == 0) {
    if (!this->squared_) {
      return p.d_grad_solP[var_][deriv_] * this->scale_;
    } else {
      return this->scale_ * 2.0 * p.grad_solP[var_][deriv_] *
             p.d_grad_solP[var_][deriv_];
    }
  } else {
    const DataType r = p.x[1];
    if (!this->squared_) {
      const DataType inv_r_comp[3] = {1. / r, 1., 1.};
      return inv_r_comp[deriv_] * p.d_grad_solP[var_][deriv_] * this->scale_;
    } else {
      const DataType rr = r * r;
      const DataType inv_rr_comp[3] = {1. / rr, 1., 1.};
      return this->scale_ * 2.0 * inv_rr_comp[deriv_] *
             p.grad_solP[var_][deriv_] * p.d_grad_solP[var_][deriv_];
    }
  }
}

} // namespace hiflow

#endif
