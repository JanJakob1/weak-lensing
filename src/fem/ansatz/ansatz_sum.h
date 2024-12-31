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

#ifndef __FEM_ANSATZ_SPACE_SUM_H_
#define __FEM_ANSATZ_SPACE_SUM_H_

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include <cmath>

#include "common/macros.h"
#include "common/vector_algebra.h"
#include "common/log.h"
#include "fem/ansatz/ansatz_space.h"


namespace hiflow {
namespace doffem {

///
/// \class AnsatzSpaceSum ansatz_sum.h
/// \brief Ansatz space that is sum of two Ansatz spaces on same reference cell and same number of components
/// \author Philipp Gerstner


template < class DataType, int DIM > 
class AnsatzSpaceSum : public virtual AnsatzSpace<DataType, DIM> 
{

public:
  typedef Vec<DIM, DataType> Coord;
  typedef boost::function< void (const Coord &, std::vector<DataType> &) > BasisEvalFunction;

  /// Default Constructor
  AnsatzSpaceSum(ConstRefCellPtr<DataType, DIM> ref_cell);

  /// Default Destructor
  virtual ~AnsatzSpaceSum();

  void init (AnsatzSpace<DataType, DIM> const * space_1, AnsatzSpace<DataType, DIM> const * space_2, AnsatzSpaceType type);  

  /// For given point, get values of all shapefunctions on reference cell
  /// The following routines implement basis functions evaluations for tensor product spaces in case nb_comp > 1.
  /// In this case, derived classes only need to provide functin evaluations for each specific component.
  /// Use the routine iv2ind for accessing function evaluations for basis function i and component var 
  // TODO: check, whether tere might be a problem with non-virtual routine in base class
  size_t iv2ind (size_t i, size_t var ) const;

  void N(const Coord &pt, std::vector< DataType > &weight) const;
  void N_x(const Coord &pt, std::vector< DataType > &weight) const;
  void N_y(const Coord &pt, std::vector< DataType > &weight) const;
  void N_z(const Coord &pt, std::vector< DataType > &weight) const;

  void N_xx(const Coord &pt, std::vector< DataType > &weight) const;
  void N_xy(const Coord &pt, std::vector< DataType > &weight) const;
  void N_xz(const Coord &pt, std::vector< DataType > &weight) const;
  void N_yy(const Coord &pt, std::vector< DataType > &weight) const;
  void N_yz(const Coord &pt, std::vector< DataType > &weight) const;
  void N_zz(const Coord &pt, std::vector< DataType > &weight) const;

protected:

  void compute_degree_hash () const 
  {
    NOT_YET_IMPLEMENTED;
  }

  void evaluate(BasisEvalFunction fun1, BasisEvalFunction fun2, const Coord &pt, std::vector< DataType > &weight) const;

  AnsatzSpace<DataType, DIM> const * space_1_;
  AnsatzSpace<DataType, DIM> const * space_2_;

  size_t dim1_;
  size_t dim2_;
  size_t weight_offset_;
};

template < class DataType, int DIM >
AnsatzSpaceSum< DataType, DIM >::AnsatzSpaceSum(ConstRefCellPtr<DataType, DIM> ref_cell)
    : AnsatzSpace<DataType, DIM>(ref_cell),
    dim1_(0),
    dim2_(0),
    space_1_(nullptr),
    space_2_(nullptr)
{
  this->name_ = "Ansatz_Sum";
  this->type_ = ANSATZ_SUM;
}

template < class DataType, int DIM > 
AnsatzSpaceSum< DataType, DIM >::~AnsatzSpaceSum() 
{
}

template < class DataType, int DIM >
void AnsatzSpaceSum< DataType, DIM >::init (AnsatzSpace<DataType, DIM> const * space_1,
                                            AnsatzSpace<DataType, DIM> const * space_2,
                                            AnsatzSpaceType type)
{
  assert (type != ANSATZ_NOT_SET);
  assert (space_1 != nullptr);
  assert (space_2 != nullptr);
  assert (space_1->dim() > 0);
  assert (space_2->dim() > 0);
  assert (space_1->nb_comp() == space_2->nb_comp());
  assert (space_1->ref_cell_type() == space_2->ref_cell_type());

  this->space_1_ = space_1;
  this->space_2_ = space_2;

  this->dim1_ = space_1->dim();
  this->dim2_ = space_2->dim();
  this->dim_ = this->dim1_ + this->dim2_;
  this->weight_offset_ = space_1->weight_size();

  this->weight_size_ = space_1->weight_size() + space_2->weight_size();
  this->nb_comp_ = space_1->nb_comp();
  this->tdim_ = space_1->tdim();

  this->max_deg_ = std::max(space_1->max_deg(), space_2->max_deg());
  this->name_ = space_1->name() + "_" + space_2->name();
  this->ref_cell_ = space_1->ref_cell();
  this->type_ = type;

  //this->compute_degree_hash();
}

template < class DataType, int DIM >
size_t AnsatzSpaceSum< DataType, DIM >::iv2ind(size_t i, size_t var) const 
{
  assert (i < this->dim1_ + this->dim2_);
  assert (this->space_1_ != nullptr);
  assert (this->space_2_ != nullptr);
  //std::cout << " sansatz sum iv " << std::endl;
  if (i < this->dim1_)
  {
//std::cout << "AnsatzSpaceSum, part 1 " << i << " , var " << var << " , index " <<  space_1_->iv2ind(i, var) << std::endl;
    return space_1_->iv2ind(i, var);
  }
  else
  { 
//std::cout << "AnsatzSpaceSum, part 2 " << i << " , var " << var << " , index " 
//              <<  this->weight_offset_ + this->space_2_->iv2ind(i-this->dim1_, var) << std::endl;
    return this->weight_offset_ + this->space_2_->iv2ind(i-this->dim1_, var);
  }
}

// TODO: avoid copy of weights by keeping N(pt, weight) from AnsatzSpace and overriding N(pt, comp, offset, weight) instead
template < class DataType, int DIM >
void AnsatzSpaceSum< DataType, DIM >::evaluate(BasisEvalFunction fun1, 
                                               BasisEvalFunction fun2, 
                                               const Coord &pt, 
                                               std::vector< DataType > &weight) const 
{ 
  assert (this->space_1_->weight_size() + this->space_2_->weight_size() == weight.size());

  std::vector<DataType> weight1 (this->space_1_->weight_size(), 0.);
  std::vector<DataType> weight2 (this->space_2_->weight_size(), 0.);

  fun1(pt, weight1);
  fun2(pt, weight2);

  const size_t size1 = weight1.size();
  const size_t size2 = weight2.size();
  
  for (size_t i = 0; i < size1; ++i)
  {
    weight[i] = weight1[i];
  }
  for (size_t i = 0; i < size2; ++i)
  {
    weight[size1+i] = weight2[i];
  }
  /* this is somehow buggy
  weight.insert(weight.begin(), weight1.begin(), weight1.end());
  weight.insert(weight.begin()+this->space_1_->weight_size(), weight2.begin(), weight2.end());
  */ 
}

template < class DataType, int DIM >
void AnsatzSpaceSum< DataType, DIM >::N(const Coord &pt, std::vector< DataType > &weight) const 
{
  BasisEvalFunction fun1 = boost::bind ( &AnsatzSpace< DataType, DIM >::N, this->space_1_, _1, _2);
  BasisEvalFunction fun2 = boost::bind ( &AnsatzSpace< DataType, DIM >::N, this->space_2_, _1, _2);
  this->evaluate(fun1, fun2, pt, weight);
}

template < class DataType, int DIM >
void AnsatzSpaceSum< DataType, DIM >::N_x(const Coord &pt, std::vector< DataType > &weight) const 
{
  BasisEvalFunction fun1 = boost::bind ( &AnsatzSpace< DataType, DIM >::N_x, this->space_1_, _1, _2);
  BasisEvalFunction fun2 = boost::bind ( &AnsatzSpace< DataType, DIM >::N_x, this->space_2_, _1, _2);
  this->evaluate(fun1, fun2, pt, weight);
}

template < class DataType, int DIM >
void AnsatzSpaceSum< DataType, DIM >::N_y(const Coord &pt, std::vector< DataType > &weight) const 
{
  BasisEvalFunction fun1 = boost::bind ( &AnsatzSpace< DataType, DIM >::N_y, this->space_1_, _1, _2);
  BasisEvalFunction fun2 = boost::bind ( &AnsatzSpace< DataType, DIM >::N_y, this->space_2_, _1, _2);
  this->evaluate(fun1, fun2, pt, weight);
}

template < class DataType, int DIM >
void AnsatzSpaceSum< DataType, DIM >::N_z(const Coord &pt, std::vector< DataType > &weight) const 
{
  BasisEvalFunction fun1 = boost::bind ( &AnsatzSpace< DataType, DIM >::N_z, this->space_1_, _1, _2);
  BasisEvalFunction fun2 = boost::bind ( &AnsatzSpace< DataType, DIM >::N_z, this->space_2_, _1, _2);
  this->evaluate(fun1, fun2, pt, weight);
}

template < class DataType, int DIM >
void AnsatzSpaceSum< DataType, DIM >::N_xx(const Coord &pt, std::vector< DataType > &weight) const 
{
  BasisEvalFunction fun1 = boost::bind ( &AnsatzSpace< DataType, DIM >::N_xx, this->space_1_, _1, _2);
  BasisEvalFunction fun2 = boost::bind ( &AnsatzSpace< DataType, DIM >::N_xx, this->space_2_, _1, _2);
  this->evaluate(fun1, fun2, pt, weight);
}

template < class DataType, int DIM >
void AnsatzSpaceSum< DataType, DIM >::N_xy(const Coord &pt, std::vector< DataType > &weight) const 
{
  BasisEvalFunction fun1 = boost::bind ( &AnsatzSpace< DataType, DIM >::N_xy, this->space_1_, _1, _2);
  BasisEvalFunction fun2 = boost::bind ( &AnsatzSpace< DataType, DIM >::N_xy, this->space_2_, _1, _2);
  this->evaluate(fun1, fun2, pt, weight);
}

template < class DataType, int DIM >
void AnsatzSpaceSum< DataType, DIM >::N_xz(const Coord &pt, std::vector< DataType > &weight) const 
{
  BasisEvalFunction fun1 = boost::bind ( &AnsatzSpace< DataType, DIM >::N_xz, this->space_1_, _1, _2);
  BasisEvalFunction fun2 = boost::bind ( &AnsatzSpace< DataType, DIM >::N_xz, this->space_2_, _1, _2);
  this->evaluate(fun1, fun2, pt, weight);
}

template < class DataType, int DIM >
void AnsatzSpaceSum< DataType, DIM >::N_yy(const Coord &pt, std::vector< DataType > &weight) const 
{
  BasisEvalFunction fun1 = boost::bind ( &AnsatzSpace< DataType, DIM >::N_yy, this->space_1_, _1, _2);
  BasisEvalFunction fun2 = boost::bind ( &AnsatzSpace< DataType, DIM >::N_yy, this->space_2_, _1, _2);
  this->evaluate(fun1, fun2, pt, weight);
}

template < class DataType, int DIM >
void AnsatzSpaceSum< DataType, DIM >::N_yz(const Coord &pt, std::vector< DataType > &weight) const 
{
  BasisEvalFunction fun1 = boost::bind ( &AnsatzSpace< DataType, DIM >::N_yz, this->space_1_, _1, _2);
  BasisEvalFunction fun2 = boost::bind ( &AnsatzSpace< DataType, DIM >::N_yz, this->space_2_, _1, _2);
  this->evaluate(fun1, fun2, pt, weight);
}

template < class DataType, int DIM >
void AnsatzSpaceSum< DataType, DIM >::N_zz(const Coord &pt, std::vector< DataType > &weight) const 
{
  BasisEvalFunction fun1 = boost::bind ( &AnsatzSpace< DataType, DIM >::N_zz, this->space_1_, _1, _2);
  BasisEvalFunction fun2 = boost::bind ( &AnsatzSpace< DataType, DIM >::N_zz, this->space_2_, _1, _2);
  this->evaluate(fun1, fun2, pt, weight);
}

} // namespace doffem
} // namespace hiflow
#endif
