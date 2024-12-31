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

#ifndef HIFLOW_VISU_FACET_VISUALIZATION
#define HIFLOW_VISU_FACET_VISUALIZATION

/// Visualization of facet within each cell, w.r.t. required material number

#include "mesh/mesh.h"

#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <sstream>
#include <tinyxml2.h>


#include <boost/function.hpp>
#include <mesh/types.h>
#include <mpi.h>

<<<<<<< HEAD:src/visualization/facet_visualization.h
#include "common/log.h"
#include "common/bbox.h"
#include "common/grid.h"
#include "fem/cell_trafo/cell_transformation.h"
#include "space/fe_evaluation.h"
#include "space/vector_space.h"
#include "mesh/attributes.h"
#include "mesh/entity.h"
#include "mesh/iterator.h"

namespace hiflow {

template < class DataType, int DIM > class VectorSpace;

/// \brief Visualization of finite element solutions.

template < class DataType, int DIM > 
class FacetVisualization {
public:
  // Type of function for evaluation.
  typedef boost::function3<
      void,
      const mesh::Entity &,            // cell
      const std::vector< Vec<DIM,DataType> > &, // reference coordinates
      std::vector< DataType > &        // values of function at the points
      >
      ScalarEvalFunction;

  typedef Vec<DIM, DataType> Coord;

  explicit FacetVisualization(const VectorSpace< DataType, DIM > &space,
                              std::vector< int > &material_numbers);

  void visualize(const FeEvalOnCellScalar<DataType, DIM> &fun, const std::string &name);
  void visualize_cell_data(const std::vector< DataType > &cell_data,
                           const std::string &name);

  void write(const std::string &filename) const;

protected:
  bool parallel_visualization_;

  const VectorSpace< DataType, DIM > &space_;
  std::map< std::string, std::vector< DataType > > functions_;
  std::map< std::string, std::vector< DataType > > functions_cell_;
  int num_visu_points_;
  int num_visu_cells_;
  std::vector< int > material_numbers_;
};

/// \brief Writer for Pvtk files.
/// \details Write PVtk files and also the corresponding Vtk files.

template < class DataType, int DIM >
class ParallelFacetVisualization : public FacetVisualization< DataType, DIM > {
public:
  /// \brief Ctor for PVtkWriter.
  /// \param [in] mpi_comm MPI Communicator.

  explicit ParallelFacetVisualization(const VectorSpace< DataType, DIM > &space,
                                      std::vector< int > &material_numbers,
                                      const MPI_Comm &mpi_comm,
                                      const int master_rank)
      : FacetVisualization< DataType, DIM >(space, material_numbers),
        comm_(mpi_comm), master_rank_(master_rank) {}

  /// \brief Writes a parallel vtk unstructured grid.
  void write(const std::string &filename, const std::string &path = "") const;

protected:
  /// The MPI Communicator.
  MPI_Comm comm_;
  const int master_rank_;
};

template < class DataType, int DIM >
FacetVisualization< DataType, DIM >::FacetVisualization( const VectorSpace< DataType, DIM > &space, std::vector< int > &material_numbers)
    : space_(space), material_numbers_(material_numbers) 
{
=======
template <class DataType>
FacetVisualization<DataType>::FacetVisualization(
  const VectorSpace<DataType> &space, const int mat_num)
  : space_(space), mat_num_(mat_num) {
>>>>>>> master:src/space/facet_visualization.cc
  parallel_visualization_ =
    space_.mesh().has_attribute("_remote_index_", space_.mesh().tdim());
}

template<class DataType>
void FacetVisualization<DataType>::get_grid_info() {
  const mesh::Mesh &mesh = this->space_.mesh();
  const mesh::TDim tdim = mesh.tdim();
  const mesh::GDim gdim = mesh.gdim();

  const mesh::TDim tdim_f = mesh.tdim()-1;
  const mesh::GDim gdim_f = mesh.gdim()-1;

  this->num_visu_points_ = 0;
  this->num_visu_cells_ = 0;

  for (mesh::EntityIterator it = mesh.begin(tdim_f), end_it = mesh.end(tdim_f);
       it != end_it; ++it) {
    if (it->get_material_number() != this->mat_num_) {
      continue;
    }

    if (this->parallel_visualization_) {
      std::vector<int> rem_ind(0);

      for (mesh::IncidentEntityIterator it_cell = it->begin_incident(tdim),
           end_it_cell = it->end_incident(tdim);
           it_cell != end_it_cell; ++it_cell) {

        int ind;
        it_cell->get<int>("_remote_index_", &ind);
        rem_ind.push_back(ind);

      }

      if (std::count(rem_ind.begin(), rem_ind.end(), -1) == 0) {
        continue;
      }
    }

    this->num_visu_cells_ += 1;
    this->num_visu_points_ += it->num_vertices();

  }
}

template <class DataType>
void FacetVisualization<DataType>::visualize_mesh() {

  this->get_grid_info();

  const mesh::Mesh &mesh = this->space_.mesh();
  const mesh::TDim tdim = mesh.tdim();
  const mesh::GDim gdim = mesh.gdim();

  const mesh::TDim tdim_f = mesh.tdim()-1;
  const mesh::GDim gdim_f = mesh.gdim()-1;

  size_t num_visu_points = this->num_visu_points_;
  size_t num_visu_cells = this->num_visu_cells_;

  ////////// collect point coordinates ////////////////////////////////////
  // determine number of points and allocate vectors for points coordinates
  this->mapped_pts_.clear();
  this->mapped_pts_.reserve(num_visu_points * 3);

  for (mesh::EntityIterator it = mesh.begin(tdim_f), end_it = mesh.end(tdim_f);
       it != end_it; ++it) {

    if (it->get_material_number() != this->mat_num_) {
      continue;
    }

    if (this->parallel_visualization_) {
      std::vector<int> rem_ind(0);

      for (mesh::IncidentEntityIterator it_cell = it->begin_incident(tdim),
           end_it_cell = it->end_incident(tdim);
           it_cell != end_it_cell; ++it_cell) {

        int ind;
        it_cell->get<int>("_remote_index_", &ind);
        rem_ind.push_back(ind);

      }

      if (std::count(rem_ind.begin(), rem_ind.end(), -1) == 0) {
        continue;
      }
    }

    std::vector<DataType> facet_coords;
    it->get_coordinates(facet_coords);

    //cell EntityIterator
    mesh::IncidentEntityIterator it_cell = it->begin_incident(tdim);

    const doffem::CellTransformation<DataType> &cell_trans =
      this->space_.GetCellTransformation(*it_cell);

    for (size_t p = 0, p_end = it->num_vertices(); p != p_end; ++p) {
      //offset on the dimension of cell
      const int offset = (gdim) * p;

      std::vector<DataType> mapped_pt(3, 0.);

      for (int c = 0; c < gdim; ++c) {
        mapped_pt[c] = facet_coords[offset + c];
      }

      this->mapped_pts_.insert(this->mapped_pts_.end(), mapped_pt.begin(),
                               mapped_pt.end());

    }
  }

  ////////// end collect point coordinates ////////////////////////////////////

  //// collect cell(for visu) data //////////////////////////////////////////////////////
  size_t p_offset = 0, cell_offset = 0;
  this->verts_.clear();
  this->verts_.reserve(num_visu_cells * 8);

  this->cell_offsets_.clear();
  this->cell_offsets_.reserve(num_visu_cells);

  this->cell_types_.clear();
  this->cell_types_.reserve(num_visu_cells);

  // Connectivity, Offsets, and Types arrays
  static const int vtk_cell_types[] = {1, 3, 5, 9};

  for (mesh::EntityIterator it = mesh.begin(tdim_f); it != mesh.end(tdim_f);
       ++it) {

    if (it->get_material_number() != this->mat_num_) {
      continue;
    }

    if (this->parallel_visualization_) {
      std::vector<int> rem_ind(0);

      for (mesh::IncidentEntityIterator it_cell = it->begin_incident(tdim),
           end_it_cell = it->end_incident(tdim);
           it_cell != end_it_cell; ++it_cell) {

        int ind;
        it_cell->get<int>("_remote_index_", &ind);
        rem_ind.push_back(ind);
      }

      if (std::count(rem_ind.begin(), rem_ind.end(), -1) == 0) {
        continue;
      }

    }

    for (int i = 0; i != it->num_vertices(); ++i) {
      this->verts_.push_back(i + cell_offset);
    }

    cell_offset += it->num_vertices();
    this->cell_offsets_.push_back(cell_offset);
    this->cell_types_.push_back(
      vtk_cell_types[static_cast<int>(it->cell_type().tag())]);

  }

}

<<<<<<< HEAD:src/visualization/facet_visualization.h
template < class DataType, int DIM >
void FacetVisualization< DataType, DIM >::visualize(const FeEvalOnCellScalar<DataType, DIM> &fun, const std::string &name) 
{
=======
template <class DataType>
void FacetVisualization<DataType>::visualize(const EvalFunction &fun,
    const std::string &name) {
>>>>>>> master:src/space/facet_visualization.cc
  const mesh::Mesh &mesh = space_.mesh();
  const mesh::TDim tdim = mesh.tdim();
  const mesh::TDim tdim_f = mesh.tdim()-1;

  if (this->mapped_pts_.size() == 0) {
    this->visualize_mesh();
  }

  std::vector<DataType> values, cell_values;
  values.reserve(this->num_visu_points_);

  for (mesh::EntityIterator it = mesh.begin(tdim_f), end_it = mesh.end(tdim_f);
       it != end_it; ++it) {

    if (it->get_material_number() != this->mat_num_) {
      continue;
    }

    if (parallel_visualization_) {
      std::vector<int> rem_ind(0);

      for (mesh::IncidentEntityIterator it_cell = it->begin_incident(tdim),
           end_it_cell = it->end_incident(tdim);
           it_cell != end_it_cell; ++it_cell) {

        int ind;
        it_cell->get<int>("_remote_index_", &ind);
        rem_ind.push_back(ind);
      }

      if (std::count(rem_ind.begin(), rem_ind.end(), -1) == 0) {
        continue;
      }

    }

<<<<<<< HEAD:src/visualization/facet_visualization.h
    const doffem::CellTransformation< DataType, DIM > &cell_trans =
        space_.get_cell_transformation(it->index());

    for (IncidentEntityIterator f_it = it->begin_incident(tdim - 1),
                                f_it_end = it->end_incident(tdim - 1);
         f_it != f_it_end; ++f_it) 
    {
=======
    cell_values.clear();
    cell_values.resize(it->num_vertices(), 1.e32);

    std::vector<DataType> coords;
    it->get_coordinates(coords);

    mesh::IncidentEntityIterator it_cell = it->begin_incident(tdim);
    const doffem::CellTransformation<DataType> &cell_trans =
      space_.GetCellTransformation(*it_cell);

    std::vector<DataType> ref_coords;
    ref_coords.resize(coords.size());

    if (tdim == 3) {
      //get ref_coords
      for (int i = 0; i != it->num_vertices(); ++i) {
        int ind_x = i * tdim + 0;
        int ind_y = i * tdim + 1;
        int ind_z = i * tdim + 2;

        cell_trans.inverse(coords[ind_x], coords[ind_y], coords[ind_z],
                           ref_coords[ind_x], ref_coords[ind_y], ref_coords[ind_z]);
>>>>>>> master:src/space/facet_visualization.cc

      }

    } else if (tdim == 2) {
      //get ref_coords
      for (int i = 0; i != it->num_vertices(); ++i) {
        int ind_x = i * tdim + 0;
        int ind_y = i * tdim + 1;

        cell_trans.inverse(coords[ind_x], coords[ind_y],
                           ref_coords[ind_x], ref_coords[ind_y]);

<<<<<<< HEAD:src/visualization/facet_visualization.h
      for (int i = 0, i_end = coords.size() / tdim; i != i_end; ++i) 
      {
        std::vector<Coord> ref_coords(1);
        Coord phys_coord; 
        phys_coord[0] = coords[i * tdim + 0];
        phys_coord[1] = coords[i * tdim + 1];
        phys_coord[2] = coords[i * tdim + 2];

        cell_trans.inverse(phys_coord, ref_coords[0]);

        fun(*it, ref_coords, cell_values);
        values.insert(values.end(), cell_values.begin(), cell_values.end());
=======
      }

    } else if (tdim == 1) {
      //get ref_coords
      for (int i = 0; i != it->num_vertices(); ++i)	 {
        int ind_x = i * tdim + 0;

        cell_trans.inverse(coords[ind_x], ref_coords[ind_x]);
>>>>>>> master:src/space/facet_visualization.cc
      }

    }

    //evalue values on each cell, and insert to values
    fun(*it_cell, ref_coords, cell_values);
    values.insert(values.end(), cell_values.begin(), cell_values.end());

  }

  functions_.insert(std::make_pair(name, values));
}

<<<<<<< HEAD:src/visualization/facet_visualization.h
template < class DataType, int DIM >
void FacetVisualization< DataType, DIM >::visualize_cell_data(const std::vector< DataType > &cell_data, const std::string &name) 
{
=======
template <class DataType>
void FacetVisualization<DataType>::visualize_cell_data(
  const std::vector<DataType> &cell_data, const std::string &name) {
>>>>>>> master:src/space/facet_visualization.cc
  const mesh::Mesh &mesh = space_.mesh();
  const mesh::TDim tdim = mesh.tdim();
  const mesh::TDim tdim_f = mesh.tdim()-1;

  if (this->mapped_pts_.size() == 0) {
    this->visualize_mesh();
  }

  std::vector<DataType> values;

  for (mesh::EntityIterator it = mesh.begin(tdim_f), end_it = mesh.end(tdim_f);
       it != end_it; ++it) {

    if (it->get_material_number() != this->mat_num_) {
      continue;
    }

    if (parallel_visualization_) {
      std::vector<int> rem_ind(0);

      for (mesh::IncidentEntityIterator it_cell = it->begin_incident(tdim),
           end_it_cell = it->end_incident(tdim);
           it_cell != end_it_cell; ++it_cell) {

        int ind;
        it_cell->get<int>("_remote_index_", &ind);
        rem_ind.push_back(ind);
      }

      if (std::count(rem_ind.begin(), rem_ind.end(), -1) == 0) {
        continue;
      }

    }

    values.insert(values.end(), 1, cell_data[it->index()]);
  }

  functions_cell_.insert(std::make_pair(name, values));
}

<<<<<<< HEAD:src/visualization/facet_visualization.h
template < class DataType, int DIM >
void FacetVisualization< DataType, DIM >::write(const std::string &filename) const 
{
  const mesh::Mesh &mesh = space_.mesh();
=======
template <class DataType>
void FacetVisualization<DataType>::write(const std::string &filename) const {

  size_t num_grid_points = this->mapped_pts_.size() / 3;
  size_t num_grid_cells = this->cell_types_.size();

  if (num_grid_points == 0) {
    LOG_INFO("Write", "No grid points visualized");
    return;
  }

  const mesh::Mesh &mesh = this->space_.mesh();
>>>>>>> master:src/space/facet_visualization.cc
  const mesh::TDim tdim = mesh.tdim();
  const mesh::GDim gdim = mesh.gdim();

  tinyxml2::XMLDocument doc;
  tinyxml2::XMLDeclaration *decl = doc.NewDeclaration();
  doc.LinkEndChild(decl);

  tinyxml2::XMLElement *vtkFile = doc.NewElement("VTKFile");
  vtkFile->SetAttribute("type", "UnstructuredGrid");
  vtkFile->SetAttribute("version", "0.1");
  vtkFile->SetAttribute("byte_order", "LittleEndian");
  vtkFile->SetAttribute("compressor", "vtkZLibDataCompressor");
  doc.LinkEndChild(vtkFile);

  tinyxml2::XMLElement *umesh = doc.NewElement("UnstructuredGrid");
  vtkFile->LinkEndChild(umesh);

  tinyxml2::XMLElement *piece = doc.NewElement("Piece");

  // TODO: checken
  piece->SetAttribute("NumberOfPoints", static_cast<int>(num_grid_points));
  piece->SetAttribute("NumberOfCells", static_cast<int>(num_grid_cells));

  umesh->LinkEndChild(piece);

  // Scratch variables for data.
  tinyxml2::XMLElement *data_array;
  std::stringstream os;

  //// Write points ////////////////////////////////////////////////////////////
  tinyxml2::XMLElement *points = doc.NewElement("Points");
  piece->LinkEndChild(points);
  data_array = doc.NewElement("DataArray");

  // Set correct length of float in dependence of DataType
  std::ostringstream type_float;
  type_float << "Float" << sizeof(DataType) * 8;
  data_array->SetAttribute("type", type_float.str().c_str());

  data_array->SetAttribute("Name", "Array");
  // always 3 comps, since vtk doesn:t handle 2D.
  data_array->SetAttribute("NumberOfComponents", "3");

  data_array->SetAttribute("format", "ascii");

  DataType range_min = std::numeric_limits<DataType>::max();
  DataType range_max = std::numeric_limits<DataType>::min();
  int cell_type_min = std::numeric_limits<int>::max();
  int cell_type_max = std::numeric_limits<int>::min();

  size_t cell_offset_max = std::numeric_limits<size_t>::min();

  for (size_t p = 0; p < num_grid_points; ++p) {
    for (int c = 0; c < 3; ++c) {
      range_min = std::min(range_min, this->mapped_pts_[3 * p + c]);
      range_max = std::max(range_max, this->mapped_pts_[3 * p + c]);

      os << this->mapped_pts_[3 * p + c] << " ";
    }
  }

  tinyxml2::XMLText *coords;
  coords = doc.NewText(os.str().c_str());

  data_array->SetAttribute("RangeMin", range_min);
  data_array->SetAttribute("RangeMax", range_max);

  points->LinkEndChild(data_array);
  data_array->LinkEndChild(coords);
  os.str("");
  os.clear();
  //// End write points
  ////////////////////////////////////////////////////////////

  //// Write cells
  /////////////////////////////////////////////////////////////////
  tinyxml2::XMLElement *cells = doc.NewElement("Cells");
  piece->LinkEndChild(cells);

  // Connectivity, Offsets, and Types arrays
  std::ostringstream off_os, type_os;

  for (size_t c = 0; c < num_grid_cells; ++c) {
    type_os << this->cell_types_[c] << " ";

    cell_type_min = std::min(cell_type_min, this->cell_types_[c]);
    cell_type_max = std::max(cell_type_max, this->cell_types_[c]);
  }

  for (size_t v = 0; v != this->verts_.size(); ++v) {
    os << this->verts_[v] << " ";
  }

  for (size_t c = 0; c != this->cell_offsets_.size(); ++c) {
    off_os << this->cell_offsets_[c] << " ";
    cell_offset_max = std::max(cell_offset_max, this->cell_offsets_[c]);
  }

  data_array = doc.NewElement("DataArray");
  data_array->SetAttribute("type", "Int64");
  data_array->SetAttribute("Name", "connectivity");
  data_array->SetAttribute("format", "ascii");
  data_array->SetAttribute("RangeMin", 0);
  data_array->SetAttribute("RangeMax", static_cast<int>(num_grid_points));

  tinyxml2::XMLText *conns = doc.NewText(os.str().c_str());
  data_array->LinkEndChild(conns);
  cells->LinkEndChild(data_array);
  os.str("");
  os.clear();

  data_array = doc.NewElement("DataArray");
  data_array->SetAttribute("type", "Int64");
  data_array->SetAttribute("Name", "offsets");
  data_array->SetAttribute("format", "ascii");
  data_array->SetAttribute("RangeMin", 0);
  data_array->SetAttribute("RangeMax", static_cast<int>(cell_offset_max));

  tinyxml2::XMLText *offs = doc.NewText(off_os.str().c_str());
  data_array->LinkEndChild(offs);
  cells->LinkEndChild(data_array);
  off_os.str("");
  off_os.clear();

  data_array = doc.NewElement("DataArray");
  data_array->SetAttribute("type", "UInt8");
  data_array->SetAttribute("Name", "types");
  data_array->SetAttribute("format", "ascii");
  data_array->SetAttribute("RangeMin", cell_type_min);
  data_array->SetAttribute("RangeMax", cell_type_max);

  tinyxml2::XMLText *types = doc.NewText(type_os.str().c_str());
  data_array->LinkEndChild(types);
  cells->LinkEndChild(data_array);
  type_os.str("");
  type_os.clear();

  //// End Write cells
  /////////////////////////////////////////////////////////////

  //// Write point data
  ////////////////////////////////////////////////////////////
  tinyxml2::XMLElement *point_data = doc.NewElement("PointData");
  piece->LinkEndChild(point_data);

  for (typename std::map<std::string, std::vector<DataType> >::const_iterator
       it = this->functions_.begin(),
       end_it = this->functions_.end();
       it != end_it; ++it) {
    data_array = doc.NewElement("DataArray");
    data_array->SetAttribute("Name", it->first.c_str());
    data_array->SetAttribute("type", type_float.str().c_str());
    data_array->SetAttribute("format", "ascii");

    for (size_t i = 0, end_i = it->second.size(); i != end_i; ++i) {
      os << (it->second)[i] << " ";
    }

    tinyxml2::XMLText *data = doc.NewText(os.str().c_str());
    data_array->LinkEndChild(data);
    point_data->LinkEndChild(data_array);
    os.str("");
    os.clear();
  }

  tinyxml2::XMLElement *cell_data = doc.NewElement("CellData");
  piece->LinkEndChild(cell_data);

  for (typename std::map<std::string, std::vector<DataType> >::const_iterator
       it = this->functions_cell_.begin(),
       end_it = this->functions_cell_.end();
       it != end_it; ++it) {

    data_array = doc.NewElement("DataArray");
    data_array->SetAttribute("Name", it->first.c_str());
    data_array->SetAttribute("type", type_float.str().c_str());
    data_array->SetAttribute("format", "ascii");

    for (size_t i = 0, end_i = it->second.size(); i != end_i; ++i) {
      os << (it->second)[i] << " ";
    }

    tinyxml2::XMLText *data = doc.NewText(os.str().c_str());
    data_array->LinkEndChild(data);
    cell_data->LinkEndChild(data_array);
    os.str("");
    os.clear();
  }

  doc.SaveFile(filename.c_str());
}

<<<<<<< HEAD:src/visualization/facet_visualization.h
//////// ParallelFacetVisualization ////////////////////////

template < class DataType, int DIM >
void ParallelFacetVisualization< DataType, DIM >::write( const std::string &filename, const std::string &path) const 
{
=======
template class FacetVisualization<double>;
template class FacetVisualization<float>;

//////// ParallelFacetVisualization ////////////////////////

template <class DataType>
void ParallelFacetVisualization<DataType>::write(
  const std::string &filename, const std::string &path_pvtu,
  const std::string &path_pvtu2path_vtu, int num_writers) const {
>>>>>>> master:src/space/facet_visualization.cc
  // const mesh::Mesh& mesh = this->space_.mesh();
  // const mesh::TDim tdim = mesh.tdim();
  // const mesh::GDim gdim = mesh.gdim();
  // get MPI rank
  int rank = -1, num_procs = -1;
  MPI_Comm_rank(comm_, &rank);
  MPI_Comm_size(comm_, &num_procs);

  std::stringstream s;
  s << rank;

  // get the correct filename including the path
  std::istringstream filename_root_dir(filename);

  std::size_t dir = filename_root_dir.str().find_last_of('.');
  LOG_DEBUG(3, "Filename: " << filename);

  std::string filename_without_suffix = filename_root_dir.str().substr(0, dir);

  assert(!filename_without_suffix.empty());

  std::string str_src_filename =
    (path_pvtu + path_pvtu2path_vtu + filename_without_suffix + "_" +
     s.str() + ".vtu");
  LOG_DEBUG(3, "Filename without suffix: " << filename_without_suffix);
  assert(!str_src_filename.empty());

<<<<<<< HEAD:src/visualization/facet_visualization.h
  // Each process writes its vtu file
  FacetVisualization< DataType, DIM >::write(str_src_filename);
=======
  int num_write_procs = num_procs;

  if (num_writers <= 0) {
    // Each process writes its vtu file
    FacetVisualization<DataType>::write(str_src_filename);

  } else {
    num_write_procs = std::min(num_writers, num_procs);

    // send data to writing procs
    this->communicate_data(num_write_procs);

    // write data
    if (rank < num_write_procs) {
      FacetVisualization<DataType>::write(str_src_filename);
    }
  }
>>>>>>> master:src/space/facet_visualization.cc

  // Master writes pvtu file
  if (rank == master_rank_) {
    tinyxml2::XMLDocument doc;
    tinyxml2::XMLDeclaration *decl = doc.NewDeclaration();
    doc.LinkEndChild(decl);

    tinyxml2::XMLElement *vtkFile = doc.NewElement("VTKFile");
    vtkFile->SetAttribute("type", "PUnstructuredGrid");
    vtkFile->SetAttribute("version", "0.1");
    vtkFile->SetAttribute("byte_order", "LittleEndian");
    vtkFile->SetAttribute("compressor", "vtkZLibDataCompressor");
    doc.LinkEndChild(vtkFile);

    tinyxml2::XMLElement *pumesh = doc.NewElement("PUnstructuredGrid");
    // GhostLevel in PUnstructuredGrid is always 0
    pumesh->SetAttribute("GhostLevel", 0);
    vtkFile->LinkEndChild(pumesh);

    tinyxml2::XMLElement *p_point_data = doc.NewElement("PPointData");
    pumesh->LinkEndChild(p_point_data);

    // Set correct length of float in dependence of DataType
    std::ostringstream type_float;
    type_float << "Float" << sizeof(DataType) * 8;

    tinyxml2::XMLElement *p_data_array;

    for (typename std::map<std::string, std::vector<DataType> >::const_iterator
         it = this->functions_.begin(),
         end_it = this->functions_.end();
         it != end_it; ++it) {

      p_data_array = doc.NewElement("PDataArray");
      p_data_array->SetAttribute("Name", it->first.c_str());
      p_data_array->SetAttribute("type", type_float.str().c_str());
      p_data_array->SetAttribute("format", "ascii");
      p_point_data->LinkEndChild(p_data_array);
    }

    tinyxml2::XMLElement *p_cell_data = doc.NewElement("PCellData");
    pumesh->LinkEndChild(p_cell_data);
    // int tdim = mesh.tdim();

    // write cell data
    // TODO currently only Float64 is supported
    tinyxml2::XMLElement *data_array;

    for (typename std::map<std::string, std::vector<DataType> >::const_iterator
         it = this->functions_cell_.begin();
         it != this->functions_cell_.end(); ++it) {
      data_array = doc.NewElement("PDataArray");
      data_array->SetAttribute("Name", it->first.c_str());
      data_array->SetAttribute("type", type_float.str().c_str());
      data_array->SetAttribute("format", "ascii");
      p_cell_data->LinkEndChild(data_array);
      p_cell_data->SetAttribute("Scalars", it->first.c_str());
    }

    // NB: This has to be AFTER the the other elements, since
    // the same order in the vtu and pvtu file is needed!

    tinyxml2::XMLElement *p_points = doc.NewElement("PPoints");
    pumesh->LinkEndChild(p_points);

    tinyxml2::XMLElement *p_points_data_array = doc.NewElement("PDataArray");
    p_points_data_array->SetAttribute("type", type_float.str().c_str());
    p_points_data_array->SetAttribute("NumberOfComponents", "3");
    p_points->LinkEndChild(p_points_data_array);

    // get the correct filename without the path
    std::size_t pos = filename_root_dir.str().find_last_of("/\\");
    assert(!filename_root_dir.str()
           .substr(pos + 1, filename_root_dir.str().length())
           .empty());

    std::stringstream str_proc_id;

    for (int proc_id = 0; proc_id < num_write_procs; ++proc_id) {
      tinyxml2::XMLElement *piece =
        doc.NewElement("Piece"); // needs to be inside the loop!
      str_proc_id << proc_id;
      std::string source_str =
        path_pvtu2path_vtu +
        filename_root_dir.str().substr(pos + 1, dir - pos - 1) + "_" +
        str_proc_id.str() + ".vtu";
      piece->SetAttribute("Source", source_str.c_str());
      pumesh->LinkEndChild(piece);
      str_proc_id.str("");
      str_proc_id.clear();
    }

    const std::string &tmp_path_pvtu = path_pvtu;
    std::string str_filename =
      (tmp_path_pvtu + filename_without_suffix + ".pvtu");
    FILE *pFile;
    pFile = fopen(str_filename.c_str(), "w");

    if (pFile != NULL) {
      doc.SaveFile(pFile);
      fclose(pFile);

    } else {
      std::stringstream err;
      err << "Path to write the files (" << str_filename << ") does not exist!";
      LOG_ERROR(err.str());
      throw std::runtime_error(err.str());
    }
  }
}

<<<<<<< HEAD:src/visualization/facet_visualization.h
=======
template <class DataType>
void ParallelFacetVisualization<DataType>::communicate_data(
  int num_write_procs) const {
  int my_rank = -1;
  int num_proc = -1;
  MPI_Comm_rank(this->mpi_comm_, &my_rank);
  MPI_Comm_size(this->mpi_comm_, &num_proc);

  assert(num_write_procs <= num_proc);

  // determine corresponding writing process for each process
  int recv_proc = -1;

  if (my_rank >= num_write_procs) {
    recv_proc = my_rank % num_write_procs;
  }

  if (recv_proc >= 0) {
    // send data to corresponding writing process
    // send point coordinates
    std::vector<size_t> send_sizes(4);
    send_sizes[0] = this->mapped_pts_.size();
    send_sizes[1] = this->verts_.size();
    send_sizes[2] = this->cell_offsets_.size();
    send_sizes[3] = this->cell_types_.size();

    std::vector<size_t> functions_sizes;
    std::vector<size_t> functions_cell_sizes;

    for (typename std::map<std::string, std::vector<DataType> >::const_iterator
         it = this->functions_.begin(),
         end_it = this->functions_.end();
         it != end_it; ++it) {
      functions_sizes.push_back(it->second.size());
    }

    for (typename std::map<std::string, std::vector<DataType> >::const_iterator
         it = this->functions_cell_.begin(),
         end_it = this->functions_cell_.end();
         it != end_it; ++it) {
      functions_cell_sizes.push_back(it->second.size());
    }

    MPI_Send(&send_sizes[0], 4, MPI_UNSIGNED_LONG, recv_proc, 0,
             this->mpi_comm_);
    MPI_Send(&functions_sizes[0], this->functions_.size(), MPI_UNSIGNED_LONG,
             recv_proc, 1, this->mpi_comm_);
    MPI_Send(&functions_cell_sizes[0], this->functions_cell_.size(),
             MPI_UNSIGNED_LONG, recv_proc, 2, this->mpi_comm_);

    MPI_Send(&this->mapped_pts_[0], this->mapped_pts_.size(), MPI_DOUBLE,
             recv_proc, 3, this->mpi_comm_);
    MPI_Send(&this->verts_[0], this->verts_.size(), MPI_INT, recv_proc, 4,
             this->mpi_comm_);
    MPI_Send(&this->cell_offsets_[0], this->cell_offsets_.size(),
             MPI_UNSIGNED_LONG, recv_proc, 5, this->mpi_comm_);
    MPI_Send(&this->cell_types_[0], this->cell_types_.size(), MPI_INT,
             recv_proc, 6, this->mpi_comm_);

    int tag = 7;

    for (typename std::map<std::string, std::vector<DataType> >::const_iterator
         it = this->functions_.begin(),
         end_it = this->functions_.end();
         it != end_it; ++it) {
      MPI_Send(&(it->second[0]), it->second.size(), MPI_DOUBLE, recv_proc, tag,
               this->mpi_comm_);
      ++tag;
    }

    for (typename std::map<std::string, std::vector<DataType> >::const_iterator
         it = this->functions_cell_.begin(),
         end_it = this->functions_cell_.end();
         it != end_it; ++it) {
      MPI_Send(&(it->second[0]), it->second.size(), MPI_DOUBLE, recv_proc, tag,
               this->mpi_comm_);
      ++tag;
    }

  } else {
    for (int p = num_write_procs; p < num_proc; ++p) {
      if (p % num_write_procs == my_rank) {
        // proc is writing process: receive data
        std::vector<size_t> send_sizes(4);

        MPI_Status status;
        MPI_Recv(&send_sizes[0], 4, MPI_UNSIGNED_LONG, p, 0, this->mpi_comm_,
                 &status);

        std::vector<size_t> functions_sizes(this->functions_.size(), 0);
        std::vector<size_t> functions_cell_sizes(this->functions_cell_.size(),
            0);

        MPI_Recv(&functions_sizes[0], this->functions_.size(),
                 MPI_UNSIGNED_LONG, p, 1, this->mpi_comm_, &status);
        MPI_Recv(&functions_cell_sizes[0], this->functions_cell_.size(),
                 MPI_UNSIGNED_LONG, p, 2, this->mpi_comm_, &status);

        std::vector<DataType> recv_pt(send_sizes[0], 0.);
        std::vector<int> recv_verts(send_sizes[1], 0);
        std::vector<size_t> recv_offsets(send_sizes[2], 0);
        std::vector<int> recv_types(send_sizes[3], 0);

        MPI_Recv(&recv_pt[0], send_sizes[0], MPI_DOUBLE, p, 3, this->mpi_comm_,
                 &status);
        MPI_Recv(&recv_verts[0], send_sizes[1], MPI_INT, p, 4, this->mpi_comm_,
                 &status);
        MPI_Recv(&recv_offsets[0], send_sizes[2], MPI_UNSIGNED_LONG, p, 5,
                 this->mpi_comm_, &status);
        MPI_Recv(&recv_types[0], send_sizes[3], MPI_INT, p, 6, this->mpi_comm_,
                 &status);

        this->mapped_pts_.insert(this->mapped_pts_.end(), recv_pt.begin(),
                                 recv_pt.end());
        this->cell_types_.insert(this->cell_types_.end(), recv_types.begin(),
                                 recv_types.end());

        size_t old_max_vert = this->verts_.size();

        for (size_t l = 0; l < recv_verts.size(); ++l) {
          recv_verts[l] += old_max_vert;
        }

        this->verts_.insert(this->verts_.end(), recv_verts.begin(),
                            recv_verts.end());

        size_t old_offset = 0;

        if (this->cell_offsets_.size() > 0) {
          old_offset = this->cell_offsets_[this->cell_offsets_.size() - 1];
        }

        for (size_t l = 0; l < recv_offsets.size(); ++l) {
          recv_offsets[l] += old_offset;
        }

        this->cell_offsets_.insert(this->cell_offsets_.end(),
                                   recv_offsets.begin(), recv_offsets.end());

        int tag = 7;
        int counter = 0;

        for (typename std::map<std::string, std::vector<DataType> >::iterator
             it = this->functions_.begin(),
             end_it = this->functions_.end();
             it != end_it; ++it) {
          size_t count = functions_sizes[counter];
          std::vector<DataType> recv_functions(count, 0.);
          MPI_Recv(&recv_functions[0], count, MPI_DOUBLE, p, tag,
                   this->mpi_comm_, &status);

          it->second.insert(it->second.end(), recv_functions.begin(),
                            recv_functions.end());
          ++counter;
          ++tag;
        }

        counter = 0;

        for (typename std::map<std::string, std::vector<DataType> >::iterator
             it = this->functions_cell_.begin(),
             end_it = this->functions_cell_.end();
             it != end_it; ++it) {
          size_t count = functions_cell_sizes[counter];
          std::vector<DataType> recv_functions_cell(count, 0.);
          MPI_Recv(&recv_functions_cell[0], count, MPI_DOUBLE, p, tag,
                   this->mpi_comm_, &status);

          it->second.insert(it->second.end(), recv_functions_cell.begin(),
                            recv_functions_cell.end());
          ++counter;
          ++tag;
        }
      }
    }
  }
}

template class ParallelFacetVisualization<double>;
template class ParallelFacetVisualization<float>;
>>>>>>> master:src/space/facet_visualization.cc
} // namespace hiflow

#endif
