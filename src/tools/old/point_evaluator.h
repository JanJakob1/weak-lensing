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

#ifndef HIFLOW_TOOLS_POINT_EVALUATOR
#define HIFLOW_TOOLS_POINT_EVALUATOR

/// \author Jonas Kratzke, Jonathan Schwegler, Simon Gawlok, Philipp Gerstner
#include "fem/cell_transformation.h"
#include "mesh/geometric_search.h"
#include "mesh/mesh.h"
#include "space/vector_space.h"
#include "space/fe_evaluation.h"
#include "mpi_tools.h"

namespace hiflow {

/// \brief Evaluates a Function in _one_ point with help of
/// the GridGeometricSearch.

template < class DataType, class T > 
class PointEvaluator
{
  // Type of function for evaluation.
  typedef boost::function3<
      void,
      const mesh::Entity &,              // cell
      const std::vector< Coordinate > &, // reference coordinates
      std::vector< T > &                 // values of function at the points
      > EvalFunction;

public:
  PointEvaluator(const VectorSpace< DataType > &space);

  /// Evaluates a Function. Returns false if the point was not found.
  /// Furthermore, the cells containing the requested point are returned.
  bool evaluate_fun(const EvalFunction &fun,
                    const std::vector< Coordinate > &point, T &value,
                    std::vector< int > &cells) const;

  bool evaluate_fun(const EvalFunction &fun,
                    const std::vector< Coordinate > &point,
                    T &value) const;
                    

  // TODO: make this function work for T instead of DataType -> need to implement appropriate mpi_allreduce routine
  
  /// Evaluates a Function and communicates the value to all processes.
  /// Returns false if point was not found (and value = 0).
  bool evaluate_fun_global(const EvalFunction &fun,
                           const std::vector< Coordinate > &point,
                           DataType &value, std::vector< int > &cells,
                           const MPI_Comm &comm) const;

  bool evaluate_fun_global(const EvalFunction &fun,
                           const std::vector< Coordinate > &point,
                           DataType &value, const MPI_Comm &comm) const;
#if 0
  void mpi_allreduce (T* send_data, T* recv_data, int count, MPI_Op op, MPI_Comm communicator ) const;
#endif

  /// Set trial cells which are searched for the given point before the complete
  /// grid is searched. Caution: Using trial cells might give a different return
  /// value if all of the following conditions are satisfied: (i) the 'point'
  /// lies on the interface of several cells (ii) not all of these cells are
  /// contained in the list of trial cells (iii) the function 'fun' is
  /// discontinuous in the given 'point'
  void set_trial_cells(const std::vector< int > &trial_cells);

private:
  const VectorSpace< DataType > &space_;
  std::vector< int > trial_cells_;
};

template < class DataType, class T >
PointEvaluator< DataType, T >::PointEvaluator(const VectorSpace< DataType > &space)
    : space_(space) 
{
  assert(space.meshPtr() != nullptr);
  trial_cells_.resize(0);
}

template < class DataType, class T >
void PointEvaluator< DataType, T >::set_trial_cells(const std::vector< int > &trial_cells) 
{
  trial_cells_.resize(trial_cells.size(), 0);
  for (size_t l = 0; l < trial_cells.size(); ++l) 
  {
    trial_cells_[l] = trial_cells[l];
  }
}

template < class DataType, class T >
bool PointEvaluator< DataType, T >::evaluate_fun( const EvalFunction &fun, 
                                                  const std::vector< Coordinate > &point,
                                                  T &value) const 
{
  std::vector< int > cell_index;
  return evaluate_fun(fun, point, value, cell_index);
}

template < class DataType, class T >
bool PointEvaluator< DataType, T >::evaluate_fun( const EvalFunction &fun, 
                                                  const std::vector< Coordinate > &point,
                                                  T &value, 
                                                  std::vector< int > &cell_index) const 
{
  T init_val {};
  value = init_val;
  mesh::MeshPtr mesh = space_.meshPtr();
  assert(mesh != 0);

  GDim gdim = mesh->gdim();

  mesh::GeometricSearchPtr search = mesh->get_geometric_search();
  assert(search != 0);

  std::vector< std::vector< Coordinate > > ref_points;

  if (!trial_cells_.empty()) 
  {
    search->find_cell(point, trial_cells_, cell_index, ref_points);
  } 
  else 
  {
    search->find_cell(point, cell_index, ref_points);
  }
  int value_count = 0;
  if (!cell_index.empty()) 
  {
    for (size_t i = 0; i < cell_index.size(); ++i) 
    {
      Entity cell = mesh->get_entity(mesh->tdim(), cell_index[i]);
      int rem_ind = -100;
      cell.get< int >("_remote_index_", &rem_ind);

      if (ref_points.size() > i && ref_points[i].size() == gdim) 
      {
        // find_cell provides reference coordinates for given cell
        std::vector< T > temp_value(1);
        fun(cell, ref_points[i], temp_value);
        value += temp_value[0];
        ++value_count;
      } 
      else 
      {
        // no reference coordinates are provided -> need to compute them
        const typename doffem::CellTransformation< DataType > &ct = space_.GetCellTransformation(cell);
        std::vector< Coordinate > ref_coord;
        bool found = ct.contains_physical_point(point, &ref_coord);
        if (found) 
        {
          std::vector< T > temp_value(1);
          fun(cell, ref_coord, temp_value);
          value += temp_value[0];
          ++value_count;
        }
      }
    }
    if (value_count > 0) 
    {
      value *= 1. / static_cast< DataType >(value_count);
      return true;
    }
    return false;
  }
  return false;
}

template < class DataType, class T >
bool PointEvaluator< DataType, T >::evaluate_fun_global(const EvalFunction &fun, 
                                                        const std::vector< Coordinate > &point,
                                                        DataType &value, 
                                                        const MPI_Comm &comm) const 
{
  std::vector< int > cell_index;
  return evaluate_fun_global(fun, point, value, cell_index, comm);
}

template < class DataType, class T >
bool PointEvaluator< DataType, T >::evaluate_fun_global( const EvalFunction &fun, 
                                                         const std::vector< Coordinate > &point,
                                                         DataType &value, 
                                                         std::vector< int > &cell_index,
                                                         const MPI_Comm &comm) const 
{
  int rank;
  MPI_Comm_rank(comm, &rank);

  DataType has_point;
  has_point = static_cast< DataType >(evaluate_fun(fun, point, value, cell_index));

  DataType sum_data[2];

  sum_data[0] = 0.;
  sum_data[1] = 0.;

  DataType origin_data[2];
  // write recv/send data in array for less communication
  origin_data[0] = value;
  origin_data[1] = has_point;
  MPI_Allreduce(origin_data, sum_data, 2, mpi_data_type< DataType >::get_type(), MPI_SUM, comm);
  if (sum_data[1] > 0.0) 
  {
    origin_data[0] = sum_data[0] / sum_data[1];
    origin_data[1] = 1.;
  } 
  else 
  {
    // no process found the point
    origin_data[0] = 0.;
    origin_data[1] = 0.;
  }
  value = origin_data[0];
  return static_cast< bool >(origin_data[1]);

#if 0
  int rank;
  MPI_Comm_rank(comm, &rank);

  DataType found_point_global = 0.;
  DataType found_point = 0.;
  found_point = static_cast< DataType >(evaluate_fun(fun, point, value, cell_index));
  
  MPI_Allreduce(&found_point, &found_point_global, 1, mpi_data_type< DataType >::get_type(), MPI_SUM, comm);
  
  T sum_data {};
  T origin_data {};
  origin_data = value;
  
  this->mpi_allreduce (&origin_data, &sum_data, 1, MPI_SUM, comm);
//  MPI_Allreduce_T <T> (&origin_data, &sum_data, 1, mpi_data_type< DataType >::get_type(), MPI_SUM, comm);
  
  if (found_point_global > 0.0) 
  {
    origin_data = sum_data * (1. / found_point_global);
    found_point = 1.;
  } 
  else 
  {
    // no process found the point
    found_point = 0.;
    origin_data *= 0.;
  }
  value = origin_data;
  return static_cast< bool >(found_point);
#endif
}

} // namespace hiflow

#endif
