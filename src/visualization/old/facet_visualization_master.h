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

#ifndef HIFLOW_SPACE_FACET_VISUALIZATION
#define HIFLOW_SPACE_FACET_VISUALIZATION

/// \author Staffan Ronnas, Martin Baumann, Teresa Beck, Simon Gawlok, Jonas
/// Kratzke
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
/// Please construct a new instantiation of the FacetVisualization every single
/// time you want to visualize your data.
///

#include "mesh/mesh.h"

#include <map>
#include <string>
#include <vector>

#include <boost/function.hpp>
#include <mesh/types.h>
#include <mpi.h>

#include "common/bbox.h"
#include "common/grid.h"
#include "space/fe_evaluation.h"

namespace hiflow {

template < class DataType > class VectorSpace;

/// \brief Visualization of finite element solutions.

template < class DataType > class FacetVisualization {
public:
  // Type of function for evaluation.
  typedef boost::function3<
  void,
  const mesh::Entity &,            // cell
  const std::vector< DataType > &, // reference coordinates
  std::vector< DataType > &        // values of function at the points
  >
  EvalFunction;

  explicit FacetVisualization(const VectorSpace< DataType > &space,
                              const int mat_num);

  ~FacetVisualization () {
    this->clear();
  }

  /// setup data structures defining the mesh in vtk format
  void visualize_mesh ();

  void visualize(const EvalFunction &fun, const std::string &name);

  void visualize_cell_data(const std::vector< DataType > &cell_data,
                           const std::string &name);

  void write(const std::string &filename) const;

  void clear() {
    this->clear_mesh_data();
    this->clear_function_data();
    this->clear_cell_data();
  }

  void clear_mesh_data() {
    this->mapped_pts_.clear();
    this->verts_.clear();
    this->cell_offsets_.clear();
    this->cell_types_.clear();
  }

  void clear_function_data() {
    this->functions_.clear();
  }

  void clear_cell_data() {
    this->functions_cell_.clear();
  }

protected:

  void get_grid_info();

  bool parallel_visualization_;

  const VectorSpace< DataType > &space_;

  int mat_num_;
  int num_visu_points_;
  int num_visu_cells_;

  mutable std::map< std::string, std::vector< DataType > > functions_;
  mutable std::map< std::string, std::vector< DataType > > functions_cell_;

  mutable std::vector<DataType> mapped_pts_;
  mutable std::vector<int> verts_;
  mutable std::vector<size_t> cell_offsets_;
  mutable std::vector<int> cell_types_;
};

/// \brief Writer for Pvtk files.
/// \details Write PVtk files and also the corresponding Vtk files.

template < class DataType >
class ParallelFacetVisualization : public FacetVisualization< DataType > {
public:
  /// \brief Ctor for PVtkWriter.
  /// \param [in] mpi_comm MPI Communicator.

  explicit ParallelFacetVisualization(const VectorSpace< DataType > &space,
                                      const int mat_num,
                                      const MPI_Comm &mpi_comm,
                                      const int master_rank)
    : FacetVisualization< DataType >(space, mat_num),
      mpi_comm_(mpi_comm), master_rank_(master_rank) {}

  /// \brief Writes a parallel vtk unstructured grid.
  /// \brief Writes a parallel vtk unstructured grid.
  void write(const std::string &filename ) const {
    this->write(filename, "", "", -1);
  }

  void write(const std::string &filename, const std::string &path_pvtu,
             const std::string &path_pvtu2path_vtu) const
  {
    this->write(filename, path_pvtu, path_pvtu2path_vtu, -1);
  }

  /// \brief Writes a parallel vtk unstructured grid.
  void write(const std::string &filename, int num_writers)  const
  {
    this->write(filename, "", "", num_writers);
  }

  /// \brief Writes a parallel vtk unstructured grid.
  void write(const std::string &filename, const std::string &path_pvtu,
             const std::string &path_pvtu2path_vtu, int num_writers) const;

protected:
  void communicate_data(int num_writers) const;

  /// The MPI Communicator.
  MPI_Comm mpi_comm_;
  const int master_rank_;



};
} // namespace hiflow

#endif
