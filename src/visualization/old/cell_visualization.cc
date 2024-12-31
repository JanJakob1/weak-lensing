#include "visualization/cell_visualization.h"

#include <cmath>
#include <limits>
#include <string>
#include <sstream>
#include <mpi.h>

#include "fem/cell_trafo/cell_transformation.h"
#include "mesh/types.h"
#include "mesh/mesh.h"
#include "mesh/attributes.h"
#include "mesh/entity.h"
#include "mesh/iterator.h"
#include "space/vector_space.h"

#include "space/element.h"

namespace hiflow {

using namespace mesh;

template <class DataType, int DIM>
CellVisualizationGrids<DataType, DIM>::CellVisualizationGrids( const Mesh *mesh,  
                                                               const int num_intervals, 
                                                               DataType origin, 
                                                               DataType side_length)
: tdim_(mesh->tdim()) 
{
  num_visu_points_ = 0;
  num_visu_cells_ = 0;
  for (mesh::EntityIterator it = mesh->begin(tdim_), end_it = mesh->end(tdim_);
       it != end_it; ++it) 
  {
    if (mesh->has_attribute("_remote_index_", tdim_)) 
    {
      int rem_ind;
      it->get<int>("_remote_index_", &rem_ind);
      if (rem_ind != -1) 
      {
#ifdef SKIP_GHOST
        continue;
#endif
      }
    }

    const CellType::Tag cell_tag = it->cell_type().tag();
    if (!grids_[cell_tag]) 
    {
      std::vector<DataType> extents(2 * tdim_);
      for (size_t i = 0; i < static_cast<size_t>(tdim_); ++i) 
      {
        extents[2 * i] = origin;
        extents[2 * i + 1] = origin + side_length;
      }
      BBox<DataType, DIM> bbox(extents);
      std::vector<int> num_intervals_vec(tdim_, num_intervals);
      grids_[cell_tag].reset(
          new Grid<DataType, DIM>(cell_tag, num_intervals_vec, bbox));
    }
    num_visu_points_ += grids_[cell_tag]->get_num_points();
    num_visu_cells_ += grids_[cell_tag]->get_num_cells();
  }
}

template <class DataType, int DIM>
int CellVisualizationGrids<DataType, DIM>::num_visu_points() const {
  return num_visu_points_;
}

template <class DataType, int DIM>
int CellVisualizationGrids<DataType, DIM>::num_visu_cells() const {
  return num_visu_cells_;
}

template <class DataType, int DIM>
int CellVisualizationGrids<DataType, DIM>::num_points(CellType::Tag cell_tag) const {
  return grids_[cell_tag]->get_num_points();
}

template <class DataType, int DIM>
int CellVisualizationGrids<DataType, DIM>::num_cells(CellType::Tag cell_tag) const {
  return grids_[cell_tag]->get_num_cells();
}

template <class DataType, int DIM>
const std::vector<int> &
CellVisualizationGrids<DataType, DIM>::vertices_of_cell(CellType::Tag cell_tag,
                                                   int i) const {
  return grids_[cell_tag]->vertices_of_cell(i);
}

template <class DataType, int DIM>
const std::vector< Vec<DIM,DataType> >&
CellVisualizationGrids<DataType, DIM>::coords(CellType::Tag cell_tag) const {
  return grids_[cell_tag]->coords();
}

template class CellVisualizationGrids<float, 1>;
template class CellVisualizationGrids<float, 2>;
template class CellVisualizationGrids<float, 3>;

template class CellVisualizationGrids<double, 1>;
template class CellVisualizationGrids<double, 2>;
template class CellVisualizationGrids<double, 3>;

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

template <class DataType, int DIM>
CellVisualization<DataType, DIM>::CellVisualization( const VectorSpace<DataType, DIM> &space, const int num_intervals)
    : Visualization<DataType, DIM> (space.get_mpi_comm()),
    space_(space), 
    grids_(&space.mesh(), num_intervals, 0., 1.) 
{
  parallel_visualization_ = space_.mesh().has_attribute("_remote_index_", space_.mesh().tdim());
  this->tdim_ = space_.mesh().tdim();
  this->gdim_ = space_.mesh().gdim();
}

template <class DataType, int DIM>
void CellVisualization<DataType, DIM>::visualize_mesh(CoordTrafoFunction const * trafo) 
{
  const mesh::Mesh &mesh = this->space_.mesh();
  const mesh::TDim tdim = mesh.tdim();
  const mesh::GDim gdim = mesh.gdim();

  size_t num_visu_points = this->grids_.num_visu_points();
  size_t num_visu_cells = this->grids_.num_visu_cells();

  ////////// collect point coordinates ////////////////////////////////////
  // determine number of points and allocate vectors for points coordinates
  this->mapped_pts_.clear();
  this->mapped_pts_.reserve(num_visu_points * 3);

  for (mesh::EntityIterator it = mesh.begin(tdim), end_it = mesh.end(tdim);
       it != end_it; ++it) 
  {
    if (this->parallel_visualization_) 
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
    const doffem::CellTransformation<DataType, DIM> &cell_trans = this->space_.get_cell_transformation(it->index());

    for (size_t p = 0, p_end = this->grids_.num_points(it->cell_type().tag());
         p != p_end; ++p) 
    {
      Coord ref_pt = this->grids_.coords(it->cell_type().tag())[p];
      Coord mapped_pt;
      cell_trans.transform(ref_pt, mapped_pt);

      if (trafo != nullptr) {
        (*trafo)(mapped_pt);
      }

      for (int d=0; d<DIM; ++d)
      {
        this->mapped_pts_.push_back(mapped_pt[d]);
      }
      for (int d=DIM; d<3; ++d)
      {
        this->mapped_pts_.push_back(0.);
      }
    }
  }

  ////////// end collect point coordinates ////////////////////////////////////

  //// collect cell data //////////////////////////////////////////////////////
  size_t p_offset = 0, cell_offset = 0;
  this->verts_.clear();
  this->verts_.reserve(num_visu_cells * 8);

  this->cell_offsets_.clear();
  this->cell_offsets_.reserve(num_visu_cells);

  this->cell_types_.clear();
  this->cell_types_.reserve(num_visu_cells);

  // Connectivity, Offsets, and Types arrays
  static const int vtk_cell_types[] = {1, 3, 5, 9, 10, 12, 14};

  for (mesh::EntityIterator it = mesh.begin(tdim); it != mesh.end(tdim); ++it) 
  {
    if (this->parallel_visualization_) 
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

    for (size_t c = 0, c_end = this->grids_.num_cells(it->cell_type().tag());
         c != c_end; ++c) 
    {
      const std::vector<int> &verts =
          this->grids_.vertices_of_cell(it->cell_type().tag(), c);
      for (size_t v = 0, v_end = verts.size(); v != v_end; ++v) 
      {
        this->verts_.push_back(verts[v] + p_offset);
      }

      cell_offset += verts.size();
      this->cell_offsets_.push_back(cell_offset);

      // pyr do not refine only into pyrs
      if (static_cast<int>(it->cell_type().tag()) == 6) 
      {
        if (c_end == 1) 
        {
          this->cell_types_.push_back(
              vtk_cell_types[static_cast<int>(it->cell_type().tag())]);
        } 
        else 
        {
          if (c < 6) 
          {
            this->cell_types_.push_back(vtk_cell_types[static_cast<int>(it->cell_type().tag())]);
          } 
          else 
          {
            this->cell_types_.push_back(vtk_cell_types[static_cast<int>(4)]);
          }
        }
      } 
      else 
      {
        this->cell_types_.push_back(vtk_cell_types[static_cast<int>(it->cell_type().tag())]);
      }
    }
    p_offset += this->grids_.num_points(it->cell_type().tag());
  }
}

template <class DataType, int DIM>
void CellVisualization<DataType, DIM>::visualize_cell_data( const std::vector<DataType> &cell_data, 
                                                            const std::string &name) 
{
  const mesh::Mesh &mesh = space_.mesh();
  const mesh::TDim tdim = mesh.tdim();

  if (this->mapped_pts_.size() == 0) 
  {
    this->visualize_mesh(nullptr);
  }

  std::vector<DataType> values;

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

    values.insert(values.end(), grids_.num_cells(it->cell_type().tag()),
                  cell_data[it->index()]);
  }

  this->functions_cell_.insert(std::make_pair(name, values));
}

template class CellVisualization<float, 1>;
template class CellVisualization<float, 2>;
template class CellVisualization<float, 3>;

template class CellVisualization<double, 1>;
template class CellVisualization<double, 2>;
template class CellVisualization<double, 3>;


} // namespace hiflow
