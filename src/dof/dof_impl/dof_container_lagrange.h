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

#ifndef _DOF_DOF_CONTAINER_LAGRANGE_H_
#define _DOF_DOF_CONTAINER_LAGRANGE_H_

#include <map>
#include <vector>

#include "dof/dof_impl/dof_container.h"

namespace hiflow {
namespace doffem {

template <class DataType, int DIM> class RefCell;
template <class DataType, int DIM> class DofFunctional;
template <class DataType, int DIM> class CellTransformation;
template <class DataType, int DIM> class FunctionSpace;
template <class DataType, int DIM> class DofPointEvaluation;

/// Predefined collection of point evaluation dof functionals for Lagrange elements 
/// on given reference cell
/// \author Philipp Gerstner

template < class DataType, int DIM > 
class DofContainerLagrange : public virtual DofContainer<DataType, DIM>
{
public:
  typedef Vec<DIM, DataType> Coord;

  /// Constructor
  DofContainerLagrange(ConstRefCellPtr<DataType, DIM> ref_cell)
  : DofContainer<DataType, DIM>(ref_cell)
  {
    this->type_ = DOF_CONTAINER_LAGRANGE;
    this->name_ = "DOF_Lagrange";
  }

  /// Destructor
  virtual ~DofContainerLagrange();

  void evaluate (FunctionSpace<DataType, DIM> const * space, 
                 const std::vector< DofID > & dof_ids, 
                 std::vector< std::vector<DataType> >& dof_values ) const;

  void evaluate (RefCellFunction<DataType, DIM> const * func, 
                 const std::vector< DofID > & dof_ids, 
                 std::vector< std::vector<DataType> >& dof_values ) const;
  
  /// initialize container for given reference cell and polynomial degrees of ansatz space
  void init (size_t degree, size_t nb_comp);
  
  std::vector< Coord > get_dof_coords () const;
  
  std::vector< Coord > get_dof_coords_on_subentity (size_t tdim, int sindex) const;
  
  
  size_t get_degree () const;

private:
  void init_dof_coord_line (size_t degree, std::vector<Coord> &dof_coord) const;
  void init_dof_coord_tri (size_t degree, std::vector<Coord> &dof_coord) const;
  void init_dof_coord_quad (size_t degree, std::vector<Coord> &dof_coord) const;
  void init_dof_coord_tet (size_t degree, std::vector<Coord> &dof_coord) const;
  void init_dof_coord_hex (size_t degree, std::vector<Coord> &dof_coord) const;
  void init_dof_coord_pyr (size_t degree, std::vector<Coord> &dof_coord) const;

  size_t degree_;
  size_t nb_comp_;
};

} // namespace doffem
} // namespace hiflow
#endif
