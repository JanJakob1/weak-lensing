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

#ifndef __FEM_FEMANAGER_H_
#define __FEM_FEMANAGER_H_

#include "common/macros.h"
#include <cassert>
#include <vector>

#include "common/vector_algebra.h"
#include "dof/dof_fem_types.h"
#include "fem/cell_transformation.h"
#include "fem/feinstance.h"
#include "fem/fetype.h"
#include "fem/femanager.h"
#include "fem/alignedhexahedrontransformation.h"
#include "fem/alignedquadtransformation.h"
#include "fem/bilinearquadtransformation.h"
#include "fem/felagrange.h"
#include "fem/felagrange_hex.h"
#include "fem/felagrange_line.h"
#include "fem/felagrange_pyr.h"
#include "fem/felagrange_quad.h"
#include "fem/felagrange_tet.h"
#include "fem/felagrange_tri.h"
#include "fem/linearlinetransformation.h"
#include "fem/linearpyramidtransformation.h"
#include "fem/lineartetrahedrontransformation.h"
#include "fem/lineartriangletransformation.h"
#include "fem/trilinearhexahedrontransformation.h"
#include "mesh/mesh.h"
#include "mesh/iterator.h"
#include "mesh/periodicity_tools.h"

namespace hiflow {
namespace doffem {

///
/// \class FEManager femanager.h
/// \brief Ancestor Manager of information about the finite element ansatz
/// \author Michael Schick<br>Martin Baumann<br>Philipp Gerstner
///

template < class DataType, int DIM >
class FEManager {
public:
  typedef Vec<DIM, DataType> Coord;
  typedef typename FEType< DataType, DIM >::ContinuityLevel ContinuityLevel;
  typedef typename FEType< DataType, DIM >::FEAnsatz FEAnsatz;

  /// Constructor with parameters geometrical dimension and nb of dof variables
  explicit FEManager(int dim, int nb_dof_var);

  /// Destructor
  virtual ~FEManager();

  /// Initialize reference to given mesh and automatically storing cell
  /// transformations
  void set_mesh(const mesh::Mesh &mesh);
  
  inline bool is_initialized () const;

  /// Information if a continuous (true) or discontinuous (false) ansatz is used
  /// for a specific variable
  inline void ccc_set_cl(int dof_var, ContinuityLevel cl);

  /// Get dimension
  inline int dim() const;

  /// Get total number of variables
  inline size_t nb_var() const;

  /// Get number of independent variables (<= get_nb_var() with "=" in case of pure tensor product FE space)
  inline size_t nb_fe() const;

  inline std::vector<int> get_phys_var(int dof_var) const;

  inline size_t var_2_fe(size_t var) const;

  /// Get the internal FE component of a physical variable (!=0 only in case of vector valued FE) 
  inline size_t var_2_comp(size_t var) const;

  /// Information if a fully continuous, tangential continuous, normal continuous or discontinuous ansatz is used
  /// for a specific variable
  inline ContinuityLevel aaa_get_cl(int dof_var) const;

  /// Information about the used FEType on a given mesh cell for a specific variable
  inline RefElement< DataType, DIM > *get_fe_on_cell_for_var(int cell_index, size_t var) const;
  inline RefElement< DataType, DIM > *get_fe_on_cell(int cell_index, size_t fe_type) const;

  /// Get the Id of a Finite Element for a specific variable
  inline typename FEType< DataType, DIM >::FEAnsatz aaa_get_fe_ansatz(int var) const;
  inline typename FEType< DataType, DIM >::FEAnsatz bbb_get_fe_ansatz(int var) const;

  inline typename FEType< DataType, DIM >::FiniteElement bbb_get_fe_type(int cell_index, int var) const;
  inline typename FEType< DataType, DIM >::FiniteElement aaa_get_fe_type(int cell_index, int dof_var) const;

  /// Get cell transformation for a specific mesh cell
  inline CellTransformation< DataType, DIM > * get_cell_transformation(int cell_index) const;

  /// Status information
  void get_status() const;

  /// Get number of initialized Finite Elements
  int nfe() const { return fe_inst_.nfe(); }

  /// \brief Initialize FE Tank by given fe ansatz and parameters for all variables (all cells get same ansatz and parameters)
  virtual void init_fe_tank(const std::vector<FEAnsatz> &ansatz, const std::vector< std::vector< int > > &param);
  
protected:
    /// \brief Initialize FE Tank by given fe ansatz and parameters for a specific variable (all cells get same ansatz)
  virtual void aaa_init_fe_tank(int var, const FEAnsatz &ansatz, const std::vector< int > &param);

  /// Initialize Linear Cell Transformations by given coordinates of mesh on
  /// cells
  void init_cell_transformation();

  void set_lagrange_element (mesh::CellType::Tag cell_type, int degree, int index);
  void init_element (std::vector<DataType> coordinates, mesh::CellType::Tag cell_type, bool is_rectangular);
  void set_cell_transformation (mesh::CellType::Tag cell_type, int align_number, int index);

  /// Dimension
  size_t dim_;

  /// Number of variables
  size_t nb_var_;

  /// Number of different FiniteElements (<= nb_var_)
  size_t nb_fe_;

  std::vector<size_t> var_2_comp_;
  std::vector<size_t> var_2_fe_;
  
  std::vector< std::vector<size_t> > fe_2_var_;

  /// This vector stores for every variable if fe_tank was initialized
  std::vector< bool > initialized_;
  bool fe_tank_initialized_;

  /// Continuous (true) or discontinuous (false) ansatz for each variable
  std::vector< bool > cl_;

  /// Const pointer to given mesh
  const mesh::Mesh *mesh_;

  /// Number of cells on the mesh
  int num_entities_;

  /// FE Tank, which holds pointers to all RefElement instances for every cell 
  std::vector< std::vector<RefElement< DataType, DIM > * > > fe_tank_;

  /// Instance holder for every needed FEType
  FEInstance< DataType, DIM > fe_inst_;

  /// Cell transformations for every mesh cell
  std::vector< CellTransformation< DataType, DIM > * > cell_transformation_;
};

//------------ INLINE FUNCTIONS FOR FEMANAGER ---------------
template < class DataType, int DIM >
inline bool FEManager< DataType, DIM >::is_initialized() const 
{
  return this->fe_tank_initialized_;
}

template < class DataType, int DIM >
inline void FEManager< DataType, DIM >::set_cont(size_t fe_ind, bool flag) 
{
  interminable_assert(fe_ind >= 0 && fe_ind < nb_fe_);
  cl_[fe_ind] = cl;
}

template < class DataType, int DIM >
inline size_t FEManager< DataType, DIM >::dim() const 
{
  return this->dim_;
}

template < class DataType, int DIM >
inline size_t FEManager< DataType, DIM >::nb_var() const 
{
  return this->nb_var_;
}

template < class DataType, int DIM >
inline size_t FEManager< DataType, DIM >::nb_fe() const 
{
  return this->nb_fe_;
}

template < class DataType, int DIM >
inline bool FEManager< DataType, DIM >::is_cont(size_t fe_ind) const 
{
  interminable_assert(fe_ind >= 0 && fe_ind < nb_fe_);
  return cl_[fe_ind];
}

template < class DataType, int DIM >
inline CellTransformation< DataType, DIM > * FEManager< DataType, DIM >::get_cell_transformation(int cell_index) const 
{
  assert(cell_index >= 0 && cell_index < cell_transformation_.size());
  return cell_transformation_[cell_index];
}

template < class DataType, int DIM >
inline size_t FEManager< DataType, DIM >::var_2_comp(size_t var) const 
{
  assert (var < var2comp_.size());
  return this->var_2_comp_[var];
}

template < class DataType, int DIM >
inline size_t FEManager< DataType, DIM >::var_2_fe(size_t var) const 
{
  assert (var < var2fe_.size());
  return this->var_2_fe_[var];
}

template < class DataType, int DIM >
inline std::vector<size_t> FEManager< DataType, DIM >::get_vars(size_t fe_ind) const 
{
  assert (fe_ind < fe2var_.size());
  return this->fe_2_var_[fe_ind];
}

template < class DataType, int DIM >
inline std::vector< RefElement< DataType, DIM > const *> FEManager< DataType, DIM >::get_fe(int cell_index) const 
{
  assert (cell_index >= 0);
  assert (cell_index < this->fe_tank_.size());

  return this->fe_tank_[cell_index];
}

template < class DataType, int DIM >
inline RefElement< DataType, DIM > const * FEManager< DataType, DIM >::get_fe(int cell_index, size_t fe_ind) const 
{
  assert (cell_index >= 0);
  assert (cell_index < this->fe_tank_.size());
  assert (fe_ind < this->fe_tank_[cell_index].size());

  return this->fe_tank_[cell_index][fe_ind];
}

template < class DataType, int DIM >
inline RefElement< DataType, DIM > * FEManager< DataType, DIM >::get_fe_for_var(int cell_index, size_t var) const 
{
  assert (var < this->nb_var_);
  return this->get_fe_on_cell(cell_index, this->var_2_fe(var));
}

template < class DataType, int DIM >
inline FEType FEManager< DataType, DIM >::fe_type_for_var(int cell_index, size_t var) const 
{
  assert (var < this->nb_var_);
  return this->get_fe_vor_var(cell_index, var)->type();
}

template < class DataType, int DIM >
inline FEType FEManager< DataType, DIM >::fe_type(int cell_index, size_t fe_ind) const 
{
  assert (feind < this->nb_fe_);
  return this->get_fe(cell_index, fe_ind)->type();
}

//------------ OTHER FUNCTIONS FOR FEMANAGER ---------------

template < class DataType, int DIM >
FEManager< DataType, DIM >::FEManager()
    : dim_(0), nb_fe_(0), nb_var_(0), num_entities_(-1), mesh_(NULL), fe_tank_initialized_(false) 
{
  initialized_.clear();

  fe_tank_.clear();
  fe_inst_.clear();

  cl_.clear()
  var_2_fe_.clear();
  var_2_comp_.clear();

}

template < class DataType, int DIM >FEManager< DataType, DIM >::~FEManager() 
{
  for (size_t i = 0, e_i = cell_transformation_.size(); i != e_i; ++i) 
  {
    delete cell_transformation_[i];
  }
}

template < class DataType, int DIM >
void FEManager< DataType, DIM >::set_mesh(const mesh::Mesh &mesh) 
{
  mesh_ = &mesh;
  num_entities_ = mesh_->num_entities(dim_);
  init_cell_transformation();
  fe_tank_.resize(num_entities_);
}

/// \details By iterating through the given mesh, for each mesh cell the
/// coordinates are extracted and an instance of a transformation with correct type
/// is stored.

template < class DataType, int DIM >
void FEManager< DataType, DIM >::init_cell_transformation() 
{
  assert(mesh_ != 0);

  cell_transformation_.clear();
  cell_transformation_.resize(num_entities_);

  std::vector< mesh::MasterSlave > period = mesh_->get_period();

  for (mesh::EntityIterator it = mesh_->begin(dim_), e_it = mesh_->end(dim_);
       it != e_it; ++it) 
  {
    // Get Cell Type (line, hex, quad, tri, tet....)

    mesh::CellType::Tag cell_type = it->cell_type().tag();
    mesh::AlignNumber align_number = it->get_align_number();

    // Initialize Cell Transformations
    this->set_cell_transformation (cell_type, align_number, it->index());  
    
    std::vector<DataType> coord_vtx;
    it->get_coordinates(coord_vtx);

    if (!period.empty()) 
    {
      std::vector<DataType> tmp_coord = coord_vtx;
      coord_vtx = unperiodify(tmp_coord, 1, period);
    }

    cell_transformation_[it->index()]->reinit(coord_vtx);
  }
}

/// \details Initializing is done by a loop over all variables

template < class DataType, int DIM >
void FEManager< DataType, DIM >::init_fe_tank( const std::vector<FEType> &fe_types, 
                                               const std::vector<bool> & continuity,
                                               const std::vector< std::vector< int > > &param) 
{
  assert(param.size() == fe_types.size());
  assert(continuity.size() == fe_types.size());

  this->cl_ = continuity;
  this->nb_fe_ = fe_types.size();

  this->var_2_fe_.clear();
  this->var_2_comp_.clear();

  this->initialized_.clear();
  this->initialized_.resize(nb_fe_, false);
  
  for (size_t fe_ind = 0; fe_ind < nb_fe_; ++fe_ind) 
  {
    // initialize FE on each mesh cell
    this->init_fe_tank(fe_ind, fe_types[fe_ind], param[fe_ind]);
    
    // setup dof variable <-> physical variable mappings 
    RefElement< DataType, DIM > * example_fe = this->fe_tank_[0][fe_ind];
    
    size_t nb_comp = example_fe->nb_comp();
    for (size_t v=0; v<nb_comp; ++v) 
    {
      this->var_2_fe_.push_back(fe_ind);
      this->var_2_comp_.push_back(v);
    }
    
    this->nb_var_ += nb_comp;
  }
  this->fe_tank_initialized_ = true;
}

/// \details By iterating through the given mesh, on each cell for a given
/// variable the
///          desired finite elements are referenced by FEInstance \see
///          FEInstance. The protected variable fe_tank_ is filled with pointers
///          to the corresponding ansatz.

template < class DataType, int DIM >
void FEManager< DataType, DIM >::init_fe_tank( size_t fe_ind, FEType fe_type, const std::vector< int > &param) 
{
  assert (fe_ind >= 0);
  assert (fe_ind < nb_fe_);
  assert(!initialized_[fe_ind]);

  assert(!param.empty());









  // set degree of continuity 
  if (ansatz == FEType< DataType, DIM >::LAGRANGE) {
    this->cl_[var] = FEType< DataType, DIM >::FULL_CONT;
  }
  else if (ansatz == FEType< DataType, DIM >::DG_LAGRANGE) {
    this->cl_[var] = FEType< DataType, DIM >::DIS_CONT;
  }
  else {
    assert(0);
  }
 
  this->fe_ansatz_[var] = ansatz;

  // Initialize pointers to existing entity types
  mesh::EntityNumber line_representer = -1;
  mesh::EntityNumber tri_representer = -1;
  mesh::EntityNumber quad_representer = -1;
  mesh::EntityNumber tet_representer = -1;
  mesh::EntityNumber hex_representer = -1;
  mesh::EntityNumber pyr_representer = -1;

  const size_t offset = var * num_entities_;
  int counter = 0;
  for (mesh::EntityIterator it = mesh_->begin(dim_), e_it = mesh_->end(dim_);
       it != e_it; ++it) {
    // 1st: Get Cell Type (Tet, Hex, Quad ...)
    mesh::CellType::Tag cell_type = it->cell_type().tag();

    int degree = param[0];
    if (param.size() == num_entities_) {
    // Assume we have user given FE degrees for each cell 
      degree = param[counter];
    }

    // 2nd: Store Finite Element in FE Tank
    if (ansatz == FEType< DataType, DIM >::LAGRANGE || ansatz == FEType< DataType, DIM >::DG_LAGRANGE) 
    {
        this->set_lagrange_element (cell_type, degree, it->index() + offset);
    }
    
    switch (cell_type) 
    {
      case mesh::CellType::LINE: 
        line_representer = it->index();
        break;
      case mesh::CellType::TRIANGLE: 
        tri_representer = it->index();
        break;
      case mesh::CellType::QUADRILATERAL: 
        quad_representer = it->index();
        break;
      case mesh::CellType::TETRAHEDRON: 
        tet_representer = it->index();
        break;
      case mesh::CellType::HEXAHEDRON: 
        hex_representer = it->index();
        break;
      case mesh::CellType::PYRAMID: 
        pyr_representer = it->index();
        break;
      default:
      assert(0);
    }
    ++counter;
  }

  // At least one cell should have been found
  assert(!(line_representer == -1 && tri_representer == -1 &&
           quad_representer == -1 && tet_representer == -1 &&
           hex_representer == -1 && pyr_representer == -1));

  // Last:  setting up numbering of reference cell via transformation of an
  // arbitrary
  //        (here the last) mesh entity back to reference cell and iteration
  //        through the vertex ids
  std::vector< DataType > coordinates;

  if (line_representer != -1) 
  {
    mesh::Entity entity = mesh_->get_entity(dim_, line_representer);
    entity.get_coordinates(coordinates);
    this->init_element (coordinates, mesh::CellType::LINE, true);
  }

  if (tri_representer != -1) 
  {
    mesh::Entity entity = mesh_->get_entity(dim_, tri_representer);
    entity.get_coordinates(coordinates);
    this->init_element (coordinates, mesh::CellType::TRIANGLE, true);
  }

  if (quad_representer != -1) 
  {
    mesh::Entity entity = mesh_->get_entity(dim_, quad_representer);
    entity.get_coordinates(coordinates);
    this->init_element (coordinates, mesh::CellType::QUADRILATERAL, this->mesh_->is_rectangular());
  }

  if (tet_representer != -1) 
  {
    mesh::Entity entity = mesh_->get_entity(dim_, tet_representer);
    entity.get_coordinates(coordinates);
    this->init_element (coordinates, mesh::CellType::TETRAHEDRON, this->mesh_->is_rectangular());
  }

  if (hex_representer != -1) 
  {
    mesh::Entity entity = mesh_->get_entity(dim_, hex_representer);
    entity.get_coordinates(coordinates);
    this->init_element (coordinates, mesh::CellType::HEXAHEDRON, this->mesh_->is_rectangular());
  }

  if (pyr_representer != -1)
  {
    mesh::Entity entity = mesh_->get_entity(dim_, pyr_representer);
    entity.get_coordinates(coordinates);
    this->init_element (coordinates, mesh::CellType::PYRAMID, this->mesh_->is_rectangular());
  }
  initialized_[var] = true;
}

template < class DataType, int DIM >
void FEManager< DataType, DIM >::get_status() const 
{
  std::cout << "DIM:   " << dim_ << std::endl;
  std::cout << "nb_fe: " << nb_fe_ << std::endl;
  std::cout << "nb_var: " << nb_var_ << std::endl;

  fe_inst_.get_status();

  std::cout << "tank size: " << fe_tank_.size() << std::endl;
  for (size_t i = 0; i < fe_tank_.size(); ++i) 
  {
    for (size_t l=0; l< fe_tank_[i].size(); ++l)
    {
      std::cout << "\t" << i << ", " << l << "\t" << fe_tank_[i][l]->get_name() << std::endl;
    }
  }
}

template <class DataType, int DIM>
void FEManager< DataType, DIM >::set_cell_transformation (mesh::CellType::Tag cell_type, int align_number, int index)
{
  switch (cell_type) 
  {
    case mesh::CellType::LINE:
      cell_transformation_[index] =
          new LinearLineTransformation< DataType, DIM >(1);
      break;
    case mesh::CellType::QUADRILATERAL:
      if (align_number >= 4) 
      {
        // 4: aligned w.r.t x and y axis
        cell_transformation_[index] =
            new AlignedQuadTransformation< DataType, DIM >(2);
      } 
      else 
      {
        cell_transformation_[index] =
            new BiLinearQuadTransformation< DataType, DIM >(2);
      }
      break;
    case mesh::CellType::TRIANGLE:
      cell_transformation_[index] =
          new LinearTriangleTransformation< DataType, DIM >(2);
      break;
    case mesh::CellType::HEXAHEDRON:
      if (align_number >= 7) 
      {
      // 7: aligned w.r.t x, y and z axis
        cell_transformation_[index] =
          new AlignedHexahedronTransformation< DataType, DIM >(3);
      } 
      else 
      { 
        cell_transformation_[index] =
          new TriLinearHexahedronTransformation< DataType, DIM >(3);
      }
      break;
    case mesh::CellType::TETRAHEDRON:
      cell_transformation_[index] =
        new LinearTetrahedronTransformation< DataType, DIM >(3);
      break;
    case mesh::CellType::PYRAMID:
      cell_transformation_[index] =
        new LinearPyramidTransformation< DataType, DIM >(3);
      break;
    default:
      assert(0);
  }
}

template <class DataType, int DIM>
void FEManager< DataType, DIM >::set_lagrange_element (mesh::CellType::Tag cell_type, int degree, int index)
{
  switch (cell_type) 
  {
    case mesh::CellType::LINE:
      this->fe_tank_[index] = fe_inst_.get_lagrange_line_element(degree);
      break;
    case mesh::CellType::QUADRILATERAL:
      this->fe_tank_[index] = fe_inst_.get_lagrange_quad_element(degree);
      break;
    case mesh::CellType::TRIANGLE:
      this->fe_tank_[index] = fe_inst_.get_lagrange_tri_element(degree);
      break;
    case mesh::CellType::HEXAHEDRON:
      this->fe_tank_[index] = fe_inst_.get_lagrange_hex_element(degree);
      break;
    case mesh::CellType::TETRAHEDRON:
      this->fe_tank_[index] = fe_inst_.get_lagrange_tet_element(degree);
      break;
    case mesh::CellType::PYRAMID:
      this->fe_tank_[index] = fe_inst_.get_lagrange_pyr_element(degree);
      break;
    default:
      assert(0);
  }
}

template <class DataType, int DIM>
void FEManager< DataType, DIM >::init_element (std::vector<DataType> coordinates, 
                                               mesh::CellType::Tag cell_type, 
                                               bool is_rectangular)
{
  int num_vert;
  switch (cell_type)
  {
    case mesh::CellType::LINE:
      num_vert = 2;
    case mesh::CellType::TRIANGLE:
      num_vert = 3;
      break;
    case mesh::CellType::QUADRILATERAL:
      num_vert = 4;
      break;
    case mesh::CellType::TETRAHEDRON:
      num_vert = 4;
      break;
    case mesh::CellType::HEXAHEDRON:
      num_vert = 8;
      break;
    case mesh::CellType::PYRAMID:
      num_vert = 5;
      break;
    default:
      assert(0);
  }
  std::vector< Vec<DIM,DataType> > ref_coords(num_vert);
  std::vector< Vec<DIM,DataType> > phys_coords(num_vert);
 
  assert (coordinates.size() == num_vert * DIM);

  for (size_t i=0; i<num_vert; ++i) 
  {
    for (size_t d=0; d<DIM; ++d) 
    {
      phys_coords[i][d] = coordinates[i*DIM+d];
    }
  }

  if (cell_type == mesh::CellType::LINE)
  {
    LinearLineTransformation<DataType, DIM> llt(1);
    llt.reinit(coordinates);

    for (size_t i = 0; i < 2; ++i) 
    {
      bool found = llt.inverse(phys_coords[i], ref_coords[i]);
      assert(found);
    }
    this->fe_inst_.init_line_elements(ref_coords);
  }
  else if (cell_type == mesh::CellType::TRIANGLE)
  {
      LinearTriangleTransformation<DataType, DIM> trafo(2);
      trafo.reinit(coordinates);
      for (size_t i = 0; i < 3; ++i) 
      {
        bool found = trafo.inverse(phys_coords[i], ref_coords[i]);
        assert(found);
      }
      this->fe_inst_.init_triangle_elements(ref_coords);
  }
  else if (cell_type == mesh::CellType::QUADRILATERAL)
  {
      if (is_rectangular)
      {
        AlignedQuadTransformation<DataType, DIM> blqt(2);
        blqt.reinit(coordinates);
        for (size_t i = 0; i < 4; ++i) 
        {
          bool found = blqt.inverse(phys_coords[i], ref_coords[i]);
          assert(found);
        }
      }
      else
      {
        BiLinearQuadTransformation<DataType, DIM> blqt(2);
        blqt.reinit(coordinates);
        for (size_t i = 0; i < 4; ++i) 
        {
          bool found = blqt.inverse(phys_coords[i], ref_coords[i]);
          assert(found);
        }
      }

      this->fe_inst_.init_quadrilateral_elements(ref_coords);
  }
  else if (cell_type == mesh::CellType::TETRAHEDRON)
  {
    LinearTetrahedronTransformation<DataType, DIM> trafo(3);
    trafo.reinit(coordinates);

    for (size_t i = 0; i < 4; ++i) 
    {
      bool found = trafo.inverse(phys_coords[i], ref_coords[i]);
      assert(found);
    }
    this->fe_inst_.init_tetrahedron_elements(ref_coords);
  }
  else if (cell_type == mesh::CellType::HEXAHEDRON)
  {
    if (is_rectangular)
    {
      AlignedHexahedronTransformation<DataType, DIM> blqt(3);
      blqt.reinit(coordinates);
      for (size_t i = 0; i < 8; ++i) 
      {
        bool found = blqt.inverse(phys_coords[i], ref_coords[i]);
        assert(found);
      }
    }
    else
    {
      TriLinearHexahedronTransformation<DataType, DIM> blqt(3);
      blqt.reinit(coordinates);
      for (size_t i = 0; i < 8; ++i) 
      {
        bool found = blqt.inverse(phys_coords[i], ref_coords[i]);
        assert(found);
      }
    }
    this->fe_inst_.init_hexahedron_elements(ref_coords);
  }
  else if (cell_type == mesh::CellType::PYRAMID)
  {
    LinearPyramidTransformation<DataType, DIM> lpr(3);
    lpr.reinit(coordinates);
    for (size_t i = 0; i < 5; ++i) 
    {
      bool found = lpr.inverse(phys_coords[i], ref_coords[i]);
      assert(found);
    }
    fe_inst_.init_pyramid_elements(ref_coords);
  }
  else
  {
      assert(0);
  }
}

} // namespace doffem
} // namespace hiflow
#endif
