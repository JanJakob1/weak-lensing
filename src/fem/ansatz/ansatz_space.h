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

#ifndef __FEM_ANSATZ_SPACE_H_
#define __FEM_ANSATZ_SPACE_H_

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include <cmath>
#include <boost/bind/bind.hpp> 

#include "common/macros.h"
#include "common/vector_algebra.h"
#include "common/log.h"
#include "fem/function_space.h"


namespace hiflow {
namespace doffem {
	
using namespace boost::placeholders;
///
/// \class AnsatzSpace ansatz_space.h
/// \brief Ancestor class of different Ansatz spaces used for defining a FE
/// \author Philipp Gerstner


template < class DataType, int DIM > 
class AnsatzSpace : public virtual FunctionSpace<DataType, DIM> 
{

public:
  typedef Vec<DIM, DataType> Coord;
  typedef boost::function< void (const Coord &, size_t, size_t, std::vector<DataType> &) > CompEvalFunction;

  /// Default Constructor
  AnsatzSpace(ConstRefCellPtr<DataType, DIM> ref_cell);

  /// Default Destructor
  virtual ~AnsatzSpace() = 0;

  inline AnsatzSpaceType type() const;

  /// Operators needed to be able to create maps where AnsatzSpace is
  /// used as key. \see FEInterfacePattern::operator < (const
  /// FEInterfacePattern& test) const Comparison by protected variable my_id_
  /// and fe_deg_ .
  virtual bool operator==(const AnsatzSpace &ansatz_slave) const;

  /// Operators needed to be able to create maps where AnsatzSpace is
  /// used as key. \see FEInterfacePattern::operator < (const
  /// FEInterfacePattern& test) const Comparison by protected variable my_id_
  /// and fe_deg_ .
  virtual bool operator<(const AnsatzSpace &ansatz_slave) const;

  /// For given point, get values of all shapefunctions on reference cell
  /// The following routines implement basis functions evaluations for tensor product spaces in case nb_comp > 1.
  /// In this case, derived classes only need to provide functin evaluations for each specific component.
  /// Use the routine iv2ind for accessing function evaluations for basis function i and component var 
  // TODO: check, whether tere might be a problem with AnsatzspaceSum
  inline size_t iv2ind (size_t i, size_t var ) const;

  /// evaluate all components of all basis functions at given point
  /// Default implementation: evaluate all components and concatenate to get multi-component evaluation
  virtual void N(const Coord &pt, std::vector< DataType > &weight) const;
  
  virtual void N_x(const Coord &pt, std::vector< DataType > &weight) const;
  virtual void N_y(const Coord &pt, std::vector< DataType > &weight) const;
  virtual void N_z(const Coord &pt, std::vector< DataType > &weight) const;

  virtual void N_xx(const Coord &pt, std::vector< DataType > &weight) const;
  virtual void N_xy(const Coord &pt, std::vector< DataType > &weight) const;
  virtual void N_xz(const Coord &pt, std::vector< DataType > &weight) const;
  virtual void N_yy(const Coord &pt, std::vector< DataType > &weight) const;
  virtual void N_yz(const Coord &pt, std::vector< DataType > &weight) const;
  virtual void N_zz(const Coord &pt, std::vector< DataType > &weight) const;

  // initialize scalar valued ansatz space of given polynomial degree 
  virtual void init (size_t degree) 
  {
    std::cerr << " init( size_t ) is not compatible with given AnsatzSpace" << std::endl;
    not_implemented();
  }
  
  // initialize vector valued ansatz space of dimension nb_comp and given uniform polynomial degree 
  virtual void init (size_t degree, size_t nb_comp) 
  {
    std::cerr << " init( size_t, size_t ) is not compatible with given AnsatzSpace" << std::endl;
    not_implemented();
  }
  
  // initialize vector valued ansatz space with given polynomial degree for each component
  virtual void init (const std::vector< size_t > &degrees)
  {
    std::cerr << " init( vector<size_t> ) is not compatible with given AnsatzSpace" << std::endl;
    not_implemented();
  }
  
  // initialize vector valued ansatz space with given polynomial degree for each component and in each space direction
  virtual void init (const std::vector< std::vector<size_t> > &degrees)
  {
    std::cerr << " init( vector< vector<size_t> > ) is not compatible with given AnsatzSpace" << std::endl;
    not_implemented();
  }
    
protected:
  virtual void compute_degree_hash () const = 0;

  /// evaluate specific component of all basis functions at given point
  /// default implementation: do nothing
  virtual void N   (const Coord &pt, size_t comp, size_t offset, std::vector<DataType> &weight) const;
  virtual void N_x (const Coord &pt, size_t comp, size_t offset, std::vector<DataType> &weight) const;
  virtual void N_y (const Coord &pt, size_t comp, size_t offset, std::vector<DataType> &weight) const;
  virtual void N_z (const Coord &pt, size_t comp, size_t offset, std::vector<DataType> &weight) const;
  virtual void N_xx(const Coord &pt, size_t comp, size_t offset, std::vector<DataType> &weight) const;
  virtual void N_xy(const Coord &pt, size_t comp, size_t offset, std::vector<DataType> &weight) const;
  virtual void N_xz(const Coord &pt, size_t comp, size_t offset, std::vector<DataType> &weight) const;
  virtual void N_yy(const Coord &pt, size_t comp, size_t offset, std::vector<DataType> &weight) const;
  virtual void N_yz(const Coord &pt, size_t comp, size_t offset, std::vector<DataType> &weight) const;
  virtual void N_zz(const Coord &pt, size_t comp, size_t offset, std::vector<DataType> &weight) const;

  AnsatzSpaceType type_;

  std::vector< size_t > comp_offset_;
  
private:
  /// helper function for combining results evaluated at diffrent components
  void evaluate(CompEvalFunction fun, const Coord &pt, std::vector< DataType > &weight) const;

};

template < class DataType, int DIM >
AnsatzSpace< DataType, DIM >::AnsatzSpace(ConstRefCellPtr<DataType, DIM> ref_cell)
    : FunctionSpace<DataType, DIM>(),
    type_(ANSATZ_NOT_SET) 
{
  this->ref_cell_ = ref_cell;
}

template < class DataType, int DIM >
AnsatzSpace< DataType, DIM >::~AnsatzSpace()
{
}

//-------------- INLINE FUNCTIONS FOR FETYPE----------------------
template < class DataType, int DIM > 
AnsatzSpaceType  AnsatzSpace< DataType, DIM >::type() const 
{
  return type_;
}

template < class DataType, int DIM >
bool AnsatzSpace< DataType, DIM >::operator==(const AnsatzSpace< DataType, DIM > &fe_slave) const 
{
  if (this->type_ == ANSATZ_NOT_SET || fe_slave.type_ == ANSATZ_NOT_SET) 
  {
    std::cout << this->type_ << " " << fe_slave.type_ << std::endl;
    assert(0);
  }
  return this->type_ == fe_slave.type_ 
      && this->nb_comp_ == fe_slave.nb_comp_
      && this->dim_ == fe_slave.dim_
      && this->deg_hash_ == fe_slave.deg_hash_;
}

template < class DataType, int DIM >
bool AnsatzSpace< DataType, DIM >::operator<(const AnsatzSpace< DataType, DIM > &slave) const 
{
  if (this->type_ == ANSATZ_NOT_SET || slave.type_ == ANSATZ_NOT_SET) {
    assert(0);
  }

  if (this->type_ < slave.type_) 
  {
    return true;
  } 
  else if (this->type_ == slave.type_) 
  {
    if (this->nb_comp_ < slave.nb_comp_)
    {
      return true;
    }
    else if (this->nb_comp_ == slave.nb_comp_)
    {
       if (this->dim_ < slave.dim_)
       {
         return true;
       }
       else if (this->dim_ == slave.dim_)
       {
         if (this->deg_hash_ < slave.deg_hash_)
         {
           return true;
         }  
       }
    }
  }
  return false;
}

/// values are stores in weights in the following order:
/// assuming a vector valued space with tensor basis
/// [psi^1_1], ............ , [psi^1_n1], ........... , [   0   ], ............ , [   0    ]
/// [   0   ], ............ , [    0   ], ........... , [   :   ], ............ , [   0    ]
/// [   :   ], ............ , [    :   ], ........... , [   0   ], ............ , [   :    ]
/// [   0   ], ............ , [    0   ], ........... , [psi^k_1], ............ , [psi^k_nk]

/// leads to
/// weights = [psi^1_1, ..., psi^1_n1, 0, ......... 0, 0, ..., 0, psi^2_1, ... , psi^2_n2, 0, ..., 0, ...........]

///           |<------ component 1 ----------------->||<--------- component 2 --------------------->||< ------...

/// basis id  |1 ........... n1        n1+1 ..... dim||1 ..... n1 n1+1 ......... n1+n2 ......., dim ||1, ....    

template < class DataType, int DIM >
size_t AnsatzSpace< DataType, DIM >::iv2ind(size_t i, size_t var) const 
{
  assert (var < this->nb_comp_);
  assert (i < this->dim_);
  return var * this->dim_ + i;
}

template < class DataType, int DIM >
void AnsatzSpace< DataType, DIM >::evaluate(CompEvalFunction fun, const Coord &pt, std::vector< DataType > &weight) const 
{
  // assume vector space is tensor product 
  assert (this->dim_ * this->nb_comp_ == weight.size());

  size_t offset = 0;

  // loop over all components
  for (size_t l=0; l<this->nb_comp_; ++l)
  {
    const size_t comp_off = this->comp_offset_[l];
    const size_t comp_size = this->comp_weight_size_[l];

    // loop over all basis functions before current component
    for (size_t i = offset; i < offset + comp_off; ++i)
    {
      weight[i] = 0.;
    } 
  
    // evaluate basis functions of current component
    fun(pt, l, offset + comp_off, weight);

    // loop over all basis functions after current component
    for (size_t i = offset + comp_off + comp_size; i < offset + this->dim_; ++i)
    {
      weight[i] = 0.;
    } 
    offset += this->dim_;
  }
}

template < class DataType, int DIM >
void AnsatzSpace< DataType, DIM >::N(const Coord &pt, std::vector< DataType > &weight) const 
{
  CompEvalFunction fun = boost::bind ( &AnsatzSpace< DataType, DIM >::N, this, _1, _2, _3, _4 );
  this->evaluate(fun, pt, weight);
//std::cout << "[AnsatzSpace] " << pt[0] << ", " << weight[0] << " ..." << std::endl;
}

template < class DataType, int DIM >
void AnsatzSpace< DataType, DIM >::N_x(const Coord &pt, std::vector< DataType > &weight) const 
{
  CompEvalFunction fun = boost::bind ( &AnsatzSpace< DataType, DIM >::N_x, this, _1, _2, _3, _4 );
  this->evaluate(fun, pt, weight);
}

template < class DataType, int DIM >
void AnsatzSpace< DataType, DIM >::N_y(const Coord &pt, std::vector< DataType > &weight) const 
{
  CompEvalFunction fun = boost::bind ( &AnsatzSpace< DataType, DIM >::N_y, this, _1, _2, _3, _4 );
  this->evaluate(fun, pt, weight);
}

template < class DataType, int DIM >
void AnsatzSpace< DataType, DIM >::N_z(const Coord &pt, std::vector< DataType > &weight) const 
{
  CompEvalFunction fun = boost::bind ( &AnsatzSpace< DataType, DIM >::N_z, this, _1, _2, _3, _4 );
  this->evaluate(fun, pt, weight);
}

template < class DataType, int DIM >
void AnsatzSpace< DataType, DIM >::N_xx(const Coord &pt, std::vector< DataType > &weight) const 
{
  CompEvalFunction fun = boost::bind ( &AnsatzSpace< DataType, DIM >::N_xx, this, _1, _2, _3, _4 );
  this->evaluate(fun, pt, weight);
}

template < class DataType, int DIM >
void AnsatzSpace< DataType, DIM >::N_xy(const Coord &pt, std::vector< DataType > &weight) const 
{
  CompEvalFunction fun = boost::bind ( &AnsatzSpace< DataType, DIM >::N_xy, this, _1, _2, _3, _4 );
  this->evaluate(fun, pt, weight);
}

template < class DataType, int DIM >
void AnsatzSpace< DataType, DIM >::N_xz(const Coord &pt, std::vector< DataType > &weight) const 
{
  CompEvalFunction fun = boost::bind ( &AnsatzSpace< DataType, DIM >::N_xz, this, _1, _2, _3, _4 );
  this->evaluate(fun, pt, weight);
}

template < class DataType, int DIM >
void AnsatzSpace< DataType, DIM >::N_yy(const Coord &pt, std::vector< DataType > &weight) const 
{
  CompEvalFunction fun = boost::bind ( &AnsatzSpace< DataType, DIM >::N_yy, this, _1, _2, _3, _4 );
  this->evaluate(fun, pt, weight);
}

template < class DataType, int DIM >
void AnsatzSpace< DataType, DIM >::N_yz(const Coord &pt, std::vector< DataType > &weight) const 
{
  CompEvalFunction fun = boost::bind ( &AnsatzSpace< DataType, DIM >::N_yz, this, _1, _2, _3, _4 );
  this->evaluate(fun, pt, weight);
}

template < class DataType, int DIM >
void AnsatzSpace< DataType, DIM >::N_zz(const Coord &pt, std::vector< DataType > &weight) const 
{
  CompEvalFunction fun = boost::bind ( &AnsatzSpace< DataType, DIM >::N_zz, this, _1, _2, _3, _4 );
  this->evaluate(fun, pt, weight);
}

template < class DataType, int DIM >
void AnsatzSpace< DataType, DIM >::N(const Coord &pt, size_t comp, size_t offset, std::vector< DataType > &weight) const 
{
  assert(comp < this->nb_comp_);
  assert(offset + this->comp_weight_size_[comp] <= weight.size());

  for (size_t ind = 0; ind < this->comp_weight_size_[comp]; ++ind)
  {
    weight[offset + ind] = 0.;
  }
}

template < class DataType, int DIM >
void AnsatzSpace< DataType, DIM >::N_x(const Coord &pt, size_t comp, size_t offset, std::vector< DataType > &weight) const 
{
  assert(comp < this->nb_comp_);
  assert(offset + this->comp_weight_size_[comp] <= weight.size());

  for (size_t ind = 0; ind < this->comp_weight_size_[comp]; ++ind)
  {
    weight[offset + ind] = 0.;
  }
}

template < class DataType, int DIM >
void AnsatzSpace< DataType, DIM >::N_y(const Coord &pt, size_t comp, size_t offset, std::vector< DataType > &weight) const 
{
  assert(comp < this->nb_comp_);
  assert(offset + this->comp_weight_size_[comp] <= weight.size());

  for (size_t ind = 0; ind < this->comp_weight_size_[comp]; ++ind)
  {
    weight[offset + ind] = 0.;
  }
}

template < class DataType, int DIM >
void AnsatzSpace< DataType, DIM >::N_z(const Coord &pt, size_t comp, size_t offset, std::vector< DataType > &weight) const 
{
  assert(comp < this->nb_comp_);
  assert(offset + this->comp_weight_size_[comp] <= weight.size());

  for (size_t ind = 0; ind < this->comp_weight_size_[comp]; ++ind)
  {
    weight[offset + ind] = 0.;
  }
}

template < class DataType, int DIM >
void AnsatzSpace< DataType, DIM >::N_xx(const Coord &pt, size_t comp, size_t offset, std::vector< DataType > &weight) const 
{
  assert(comp < this->nb_comp_);
  assert(offset + this->comp_weight_size_[comp] <= weight.size());

  for (size_t ind = 0; ind < this->comp_weight_size_[comp]; ++ind)
  {
    weight[offset + ind] = 0.;
  }
}

template < class DataType, int DIM >
void AnsatzSpace< DataType, DIM >::N_xy(const Coord &pt, size_t comp, size_t offset, std::vector< DataType > &weight) const 
{
  assert(comp < this->nb_comp_);
  assert(offset + this->comp_weight_size_[comp] <= weight.size());

  for (size_t ind = 0; ind < this->comp_weight_size_[comp]; ++ind)
  {
    weight[offset + ind] = 0.;
  }
}

template < class DataType, int DIM >
void AnsatzSpace< DataType, DIM >::N_xz(const Coord &pt, size_t comp, size_t offset, std::vector< DataType > &weight) const 
{
  assert(comp < this->nb_comp_);
  assert(offset + this->comp_weight_size_[comp] <= weight.size());

  for (size_t ind = 0; ind < this->comp_weight_size_[comp]; ++ind)
  {
    weight[offset + ind] = 0.;
  }
}

template < class DataType, int DIM >
void AnsatzSpace< DataType, DIM >::N_yy(const Coord &pt, size_t comp, size_t offset, std::vector< DataType > &weight) const 
{
  assert(comp < this->nb_comp_);
  assert(offset + this->comp_weight_size_[comp] <= weight.size());

  for (size_t ind = 0; ind < this->comp_weight_size_[comp]; ++ind)
  {
    weight[offset + ind] = 0.;
  }
}

template < class DataType, int DIM >
void AnsatzSpace< DataType, DIM >::N_yz(const Coord &pt, size_t comp, size_t offset, std::vector< DataType > &weight) const 
{
  assert(comp < this->nb_comp_);
  assert(offset + this->comp_weight_size_[comp] <= weight.size());

  for (size_t ind = 0; ind < this->comp_weight_size_[comp]; ++ind)
  {
    weight[offset + ind] = 0.;
  }
}

template < class DataType, int DIM >
void AnsatzSpace< DataType, DIM >::N_zz(const Coord &pt, size_t comp, size_t offset, std::vector< DataType > &weight) const 
{
  assert(comp < this->nb_comp_);
  assert(offset + this->comp_weight_size_[comp] <= weight.size());

  for (size_t ind = 0; ind < this->comp_weight_size_[comp]; ++ind)
  {
    weight[offset + ind] = 0.;
  }
}

} // namespace doffem
} // namespace hiflow
#endif
