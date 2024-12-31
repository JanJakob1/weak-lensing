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

#ifndef HIFLOW_VISU_CELL_VISUALIZATION
#define HIFLOW_VISU_CELL_VISUALIZATION

/// \author Staffan Ronnas, Martin Baumann, Teresa Beck, Simon Gawlok, Jonas
/// Kratzke, Philipp Gerstner
///
/// \brief Visualization of finite element functions.
///
/// Using this class a Vtk (http://www.vtk.org/) unstructured grid visualization
/// file can be created. Please find detailed information about Vtk's file
/// formats at http://www.vtk.org/VTK/img/file-formats.pdf.
/// This type of visualization writes out every cell and with function values
/// provided by a user-defined evaluation function.
///
/// Please note for simulations with multiple visualization calls, that this
/// class is NOT ment to be initialized once for several visualization calls.
/// Please construct a new instantiation of the CellVisualization every single
/// time you want to visualize your data.
///


#include <map>
#include <string>
#include <vector>
#include <sstream>
#include <boost/function.hpp>

#include "visualization/visualization.h"
#include "common/log.h"
#include "common/bbox.h"
#include "common/grid.h"
#include "mesh/iterator.h"

#define nSKIP_GHOST

namespace hiflow {

using namespace mesh;

template < class DataType, int DIM > class VectorSpace;

/// \brief Description of a square 2d grid or cubic 3d grid.

template < class DataType, int DIM > 
class CellVisualizationGrids {
public:
  typedef Vec<DIM, DataType> Coord;

  CellVisualizationGrids(const Mesh *mesh, int num_intervals, DataType origin, DataType side_length);

  int num_visu_points() const;
  int num_visu_cells() const;
  int num_points(CellType::Tag cell_tag) const;
  int num_cells(CellType::Tag cell_tag) const;

  const std::vector< int > &vertices_of_cell(CellType::Tag cell_tag, int i) const;
  const std::vector< Coord > &coords(CellType::Tag cell_tag) const;

private:
  typename ScopedPtr< Grid< DataType, DIM > >::Type grids_[CellType::NUM_CELL_TYPES];
  const int tdim_;
  int num_visu_points_;
  int num_visu_cells_;
};

/// \brief Visualization of finite element solutions.

template < class DataType, int DIM > 
class CellVisualization : public virtual Visualization<DataType, DIM>
{
public:
  typedef boost::function1<void, Vec<DIM, DataType>&> CoordTrafoFunction;
  
  typedef Vec<DIM, DataType> Coord;

  explicit CellVisualization(const VectorSpace< DataType, DIM > &space, const int num_intervals);

  ~CellVisualization () {
    this->clear();
  }
  
  /// setup data structures defining the mesh in vtk format
  void visualize_mesh (CoordTrafoFunction const * trafo = nullptr);
  
  template <class CellWiseEvaluator>
  void visualize (const CellWiseEvaluator &fun, const std::vector<std::string> &name);

  template <class CellWiseEvaluator>
  void visualize_grad (const CellWiseEvaluator &fun, const std::vector<std::string> &name);

  void visualize_cell_data(const std::vector< DataType > &cell_data, const std::string &name);

protected:
  bool parallel_visualization_;
        
  const VectorSpace< DataType, DIM > &space_;
  CellVisualizationGrids< DataType, DIM > grids_;
  
};

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

template <class DataType, int DIM>
template <class CellWiseEvaluator>
void CellVisualization<DataType, DIM>::visualize(const CellWiseEvaluator &fun, 
                                                 const std::vector<std::string> &names)
{
  const mesh::Mesh &mesh = space_.mesh();
  const mesh::TDim tdim = mesh.tdim();

  if (this->mapped_pts_.size() == 0) 
  {
    this->visualize_mesh(nullptr);
  }

  const size_t nb_var = fun.nb_comp();
  assert (nb_var == names.size());
  
  std::vector< std::vector<DataType> > values(nb_var); 
  for (size_t d=0; d<nb_var; ++d)
  {
    values[d].reserve(grids_.num_visu_points());
  }
  
  std::vector< DataType > cell_values(nb_var, 0.); 
  
  for (mesh::EntityIterator it = mesh.begin(tdim), end_it = mesh.end(tdim);
       it != end_it; ++it) 
  {

    if (parallel_visualization_) 
    {
      int rem_ind = -100;
      it->get<int>("_remote_index_", &rem_ind);
      if (rem_ind != -1) 
      {
#ifdef SKIP_GHOST
        continue;
#endif
      }
    }

    std::vector< Coord > ref_pts = this->grids_.coords(it->cell_type().tag());
    const size_t num_pt = ref_pts.size();
    
    /// loop over points in cell
    for (size_t i=0; i<num_pt; ++i)
    {
      fun.r_evaluate(*it, ref_pts[i], cell_values); 
      assert (cell_values.size() == nb_var);
    
      for (size_t d=0; d<nb_var; ++d)
      {
        values[d].push_back(cell_values[d]);
      }
    }
  }
  for (size_t d=0; d<nb_var; ++d)
  {
    this->functions_.insert(std::make_pair(names[d], values[d]));
  }
}


template <class DataType, int DIM>
template <class CellWiseEvaluator>
void CellVisualization<DataType, DIM>::visualize_grad(const CellWiseEvaluator &fun, 
                                                 const std::vector<std::string> &names)
{
  const mesh::Mesh &mesh = space_.mesh();
  const mesh::TDim tdim = mesh.tdim();

  if (this->mapped_pts_.size() == 0) 
  {
    this->visualize_mesh(nullptr);
  }

  const size_t nb_var = fun.nb_comp();
  assert (nb_var == names.size());
  
  std::vector<  std::vector<Vec<DIM, DataType > > > values(nb_var); //
  
  for (size_t d=0; d<nb_var; ++d)
  {
      values[d].reserve(grids_.num_visu_points());
  }
  
  std::vector< Vec<DIM, DataType> > cell_values(nb_var); //
  
  
  for (mesh::EntityIterator it = mesh.begin(tdim), end_it = mesh.end(tdim);
       it != end_it; ++it) 
  {

    if (parallel_visualization_) 
    {
      int rem_ind = -100;
      it->get<int>("_remote_index_", &rem_ind);
      if (rem_ind != -1) 
      {
#ifdef SKIP_GHOST
        continue;
#endif
      }
    }

    std::vector< Coord > ref_pts = this->grids_.coords(it->cell_type().tag());
    const size_t num_pt = ref_pts.size();
    
    /// loop over points in cell
    for (size_t i=0; i<num_pt; ++i)
    {
      fun.r_evaluate_grad(*it, ref_pts[i], cell_values); //
      assert (cell_values.size() == nb_var);
    
      for (size_t d=0; d<nb_var; ++d)
      {
          values[d].push_back(cell_values[d]);
      }
    }
  }
  for (size_t d=0; d<nb_var; ++d)
  {
      this->functions_grad_.insert(std::make_pair(names[d], values[d]));
  }
}






} // namespace hiflow

#endif
