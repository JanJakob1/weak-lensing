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

#include "fem/fe_instance.h"
#include "fem/fe_reference.h"
#include "fem/fe_transformation.h"
#include "fem/reference_cell.h"
#include "dof/dof_impl/dof_container.h"
#include "dof/dof_impl/dof_container_rt_bdm.h"
#include "dof/dof_impl/dof_container_lagrange.h"
#include "fem/ansatz/ansatz_space.h"
#include "fem/ansatz/ansatz_sum.h"
#include "fem/ansatz/ansatz_p_line_lagrange.h"
#include "fem/ansatz/ansatz_p_tri_lagrange.h"
#include "fem/ansatz/ansatz_aug_p_tri_mono.h"
#include "fem/ansatz/ansatz_p_tet_lagrange.h"
//#include "fem/ansatz/ansatz_pyr_lagrange.h"
#include "fem/ansatz/ansatz_q_quad_lagrange.h"
#include "fem/ansatz/ansatz_q_hex_lagrange.h"
#include "fem/cell_trafo/cell_transformation.h"
#include "fem/cell_trafo/linearlinetransformation.h"
#include "fem/cell_trafo/lineartriangletransformation.h"
#include "fem/cell_trafo/lineartetrahedrontransformation.h"
#include "fem/cell_trafo/linearpyramidtransformation.h"
#include "fem/cell_trafo/bilinearquadtransformation.h"
#include "fem/cell_trafo/trilinearhexahedrontransformation.h"
#include "fem/cell_trafo/alignedquadtransformation.h"
#include "fem/cell_trafo/alignedhexahedrontransformation.h"
#include "mesh/geometric_tools.h"

#define nFORCE_NONLINEAR_TRAFO

namespace hiflow {
namespace doffem {

template <class DataType, int DIM>
void create_RTBDM_container (DofContainerType type, 
                             size_t deg, 
                             ConstRefCellPtr<DataType, DIM> ref_cell, 
                             DofContainer<DataType, DIM>** dofs);

template <>
void create_RTBDM_container (DofContainerType type, 
                             size_t deg, 
                             ConstRefCellPtr<float, 1> ref_cell, 
                             DofContainer<float, 1>** dofs) 
{
  std::cout << "RT BDM elements are not available for 1D " << std::endl;
  quit_program();
}

template <>
void create_RTBDM_container (DofContainerType type, 
                             size_t deg, 
                             ConstRefCellPtr<double, 1>  ref_cell, 
                             DofContainer<double, 1>** dofs) 
{
  std::cout << "RT BDM elements are not available for 1D " << std::endl;
  quit_program();
}

template <>
void create_RTBDM_container (DofContainerType type, 
                             size_t deg, 
                             ConstRefCellPtr<float, 2>  ref_cell, 
                             DofContainer<float, 2>** dofs) 
{
  // create DofContainer corresponding to BDM elements
  *dofs = new DofContainerRTBDM<float, 2>(ref_cell);
  dynamic_cast<DofContainerRTBDM<float, 2> *>(*dofs)->init(deg, type);
}

template <>
void create_RTBDM_container (DofContainerType type, 
                             size_t deg, 
                             ConstRefCellPtr<double, 2>  ref_cell, 
                             DofContainer<double, 2>** dofs) 
{
  // create DofContainer corresponding to BDM elements
  *dofs = new DofContainerRTBDM<double, 2>(ref_cell);
  dynamic_cast<DofContainerRTBDM<double, 2> *>(*dofs)->init(deg, type);
}

template <>
void create_RTBDM_container (DofContainerType type, 
                             size_t deg, 
                             ConstRefCellPtr<float, 3>  ref_cell, 
                             DofContainer<float, 3>** dofs) 
{
  // create DofContainer corresponding to BDM elements
  *dofs = new DofContainerRTBDM<float, 3>(ref_cell);
  dynamic_cast<DofContainerRTBDM<float, 3> *>(*dofs)->init(deg, type);
}

template <>
void create_RTBDM_container (DofContainerType type, 
                             size_t deg, 
                             ConstRefCellPtr<double, 3>  ref_cell, 
                             DofContainer<double, 3>** dofs) 
{
  // create DofContainer corresponding to BDM elements
  *dofs = new DofContainerRTBDM<double, 3>(ref_cell);
  dynamic_cast<DofContainerRTBDM<double, 3> *>(*dofs)->init(deg, type);
}

// degrees[c][d] : polynomial degree of component c in spatial direction d
template <class DataType, int DIM>
void create_lagrange_element (mesh::CellType::Tag topo_cell_type,
                              size_t degree, 
                              size_t nb_comp,
                              AnsatzSpace<DataType, DIM>** ansatz,
                              DofContainer<DataType, DIM>** dofs,
                              ConstRefCellPtr<DataType, DIM>& ref_cell,
                              FETransformation<DataType, DIM>** fe_trafo,
                              RefElement<DataType, DIM>* ref_elem,
                              bool& is_modal_basis,
                              FEConformity& conform) 
{
 
  // create ansatz space object and select referenc cell
  switch (topo_cell_type)
  {
    case mesh::CellType::LINE:
      ref_cell = ConstRefCellPtr<DataType, DIM>(new RefCellLineStd<DataType, DIM>());
      *ansatz = new PLineLag<DataType, DIM>(ref_cell);
      break;
    case mesh::CellType::TRIANGLE:
      ref_cell = ConstRefCellPtr<DataType, DIM>(new RefCellTriStd<DataType, DIM>());
      *ansatz = new PTriLag<DataType, DIM>(ref_cell);
      break;
    case mesh::CellType::QUADRILATERAL:
      ref_cell = ConstRefCellPtr<DataType, DIM>(new RefCellQuadStd<DataType, DIM>());
      *ansatz = new QQuadLag<DataType, DIM>(ref_cell);
      break;
    case mesh::CellType::TETRAHEDRON:
      ref_cell = ConstRefCellPtr<DataType, DIM>(new RefCellTetStd<DataType, DIM>());
      *ansatz = new PTetLag<DataType, DIM>(ref_cell);
      break;
    case mesh::CellType::HEXAHEDRON:
      ref_cell = ConstRefCellPtr<DataType, DIM>(new RefCellHexStd<DataType, DIM>());
      *ansatz = new QHexLag<DataType, DIM>(ref_cell);
      break;
    case mesh::CellType::PYRAMID:
      NOT_YET_IMPLEMENTED;
      ref_cell = ConstRefCellPtr<DataType, DIM>(new RefCellPyrStd<DataType, DIM>());
//    *ansatz = new PyrLag<DataType, DIM>(ref_cell);
      break;
    default:
      std::cout << "Unexpected reference cell type " << std::endl;
      quit_program();
  }
  (*ansatz)->init(degree, nb_comp);
  
  // create DofContainer correpsonding to Lagrange elements
  *dofs = new DofContainerLagrange<DataType, DIM>(ref_cell);
  dynamic_cast<DofContainerLagrange<DataType, DIM> * > (*dofs)->init(degree, nb_comp);
      
  // create FE transformation
  *fe_trafo = new FETransformationStandard<DataType, DIM>(ref_elem);
  
  is_modal_basis = true;
  conform = FE_CONFORM_H1;
}

template <class DataType, int DIM>
void create_BDM_element (mesh::CellType::Tag topo_cell_type, 
                         size_t deg, 
                         AnsatzSpace<DataType, DIM>** ansatz,
                         DofContainer<DataType, DIM>** dofs,
                         ConstRefCellPtr<DataType, DIM>& ref_cell,
                         FETransformation<DataType, DIM>** fe_trafo,
                         RefElement<DataType, DIM>* ref_elem,
                         bool& is_modal_basis,
                         FEConformity& conform) 
{
  // create ansatz space object     
  switch (topo_cell_type)
  {
    case mesh::CellType::TRIANGLE:
      ref_cell = ConstRefCellPtr<DataType, DIM>(new RefCellTriStd<DataType, DIM>());
      *ansatz = new PTriLag<DataType, DIM>(ref_cell);
      break;
    case mesh::CellType::QUADRILATERAL:
      NOT_YET_IMPLEMENTED;
      break;
    case mesh::CellType::TETRAHEDRON:
      NOT_YET_IMPLEMENTED;
      break;
    case mesh::CellType::HEXAHEDRON:
      NOT_YET_IMPLEMENTED;
      break;
    default:
      std::cout << "Unexpected reference cell type " << std::endl;
      quit_program();
  }
      
  (*ansatz)->init(deg, DIM);
  
  create_RTBDM_container (DOF_CONTAINER_BDM, deg, ref_cell, dofs);
   
  // create e transformation
  *fe_trafo = new FETransformationContraPiola<DataType, DIM>(ref_elem);
    
  is_modal_basis = false;
  conform = FE_CONFORM_HDIV;
}

template <class DataType, int DIM>
void create_RT_element (mesh::CellType::Tag topo_cell_type, 
                         size_t deg, 
                         AnsatzSpace<DataType, DIM>** ansatz,
                         DofContainer<DataType, DIM>** dofs,
                         ConstRefCellPtr<DataType, DIM>& ref_cell,
                         FETransformation<DataType, DIM>** fe_trafo,
                         RefElement<DataType, DIM>* ref_elem,
                         bool& is_modal_basis,
                         FEConformity& conform) 
{
  // create ansatz space object     
  if (topo_cell_type == mesh::CellType::TRIANGLE)
  {
    ref_cell = ConstRefCellPtr<DataType, DIM> (new RefCellTriStd<DataType, DIM>());
      
    PTriLag<DataType, DIM>* ansatz1 = new PTriLag<DataType, DIM>(ref_cell);
    AugPTriMono<DataType, DIM>* ansatz2 = new AugPTriMono<DataType, DIM>(ref_cell);
    ansatz1->init(deg, DIM);
    ansatz2->init(deg);
    AnsatzSpaceSum<DataType, DIM> * ansatz_sum = new AnsatzSpaceSum<DataType, DIM>(ref_cell);
    ansatz_sum->init(ansatz1, ansatz2, ANSATZ_RT);
    
    *ansatz = ansatz_sum;
  }
  else if (topo_cell_type == mesh::CellType::QUADRILATERAL)
  {
    NOT_YET_IMPLEMENTED;
  }
  else if (topo_cell_type == mesh::CellType::TETRAHEDRON)
  {
    NOT_YET_IMPLEMENTED;
  }
  else if (topo_cell_type == mesh::CellType::HEXAHEDRON)
  {
    NOT_YET_IMPLEMENTED;
  }
  else 
  {
    std::cout << "Unexpected reference cell type " << std::endl;
    quit_program();
  }
  
  create_RTBDM_container (DOF_CONTAINER_RT, deg, ref_cell, dofs);
   
  // create e transformation
  *fe_trafo = new FETransformationContraPiola<DataType, DIM>(ref_elem);
    
  is_modal_basis = false;
  conform = FE_CONFORM_HDIV;
}


template <class DataType, int DIM>
void FEInstance<DataType, DIM>::clear()
{
  for (size_t i=0; i<this->ref_elements_.size(); ++i)
  {
    if (this->ref_elements_[i] != nullptr)
    {
      delete this->ref_elements_[i];
      this->ref_elements_[i] = nullptr;
    }
    if (this->ansatz_spaces_[i] != nullptr)
    {
      delete this->ansatz_spaces_[i];
      this->ansatz_spaces_[i] = nullptr;
    }
    if (this->dof_containers_[i] != nullptr)
    {
      delete this->dof_containers_[i];
      this->dof_containers_[i] = nullptr;
    }
    if (this->fe_trafos_[i] != nullptr)
    {
      delete this->fe_trafos_[i];
      this->fe_trafos_[i] = nullptr;
    }
    /*
    if (this->ref_cells_[i] != nullptr)
    {
      delete this->ref_cells_[i];
      this->ref_cells_[i] = nullptr;
    }
    */
  }
  this->lagrange_only_ = true;
  this->ref_elements_.clear();
  this->ansatz_spaces_.clear();
  this->dof_containers_.clear();
  this->fe_trafos_.clear();
  this->max_fe_conform_.clear();
  this->nb_nonlin_trafos_ = 0;
}

template <class DataType, int DIM>
RefElement<DataType, DIM> const * FEInstance<DataType, DIM>::get_fe (size_t fe_id) const
{
  assert (fe_id < this->ref_elements_.size());
  return this->ref_elements_[fe_id];
}

template <class DataType, int DIM>
size_t FEInstance<DataType, DIM>::add_fe ( FEType fe_type, 
                                           mesh::CellType::Tag topo_cell_type, 
                                           const std::vector<int> &param )
{
  // create instance of a Finite Element, defined by general type, topology of reference cell and additional parameters
  
  // first check whether instance was already created. If yes, return corresponding index
  for (size_t l=0; l<added_fe_types_.size(); ++l)
  {
    if (fe_type == added_fe_types_[l])
    {
      if (topo_cell_type == added_cell_types_[l])
      {
        if (param == added_params_[l])
        {
          return l;
        }
      }
    }
  }
  
  AnsatzSpace<DataType, DIM>* ansatz = nullptr;
  DofContainer<DataType, DIM>* dofs = nullptr;
  FETransformation<DataType, DIM>* fe_trafo = nullptr;
  ConstRefCellPtr<DataType, DIM> ref_cell;
  RefElement<DataType, DIM>* ref_elem = new RefElement<DataType, DIM>();
  RefCellType cell_type;
    
  // Note: modal_basis = true <=> the dofs sigma_i and the basis functions phi_j in ansatz satisfy dof_i(phi_j) = delta_{ij} by construction
  bool modal_basis = false;
  FEConformity conform = FE_CONFORM_L2;
  
  // This is the only location that has to be modified after having implemented a new type of element
  if (fe_type == FE_TYPE_LAGRANGE)
  {
    assert (param.size() == 1);
    assert (param[0] >= 0);
    create_lagrange_element(topo_cell_type, param[0], 1, &ansatz, &dofs, ref_cell, &fe_trafo, ref_elem, modal_basis, conform);
  }
  else if (fe_type == FE_TYPE_LAGRANGE_VECTOR)
  {
    assert (param.size() == 1);
    assert (param[0] >= 0);
    create_lagrange_element(topo_cell_type, param[0], DIM, &ansatz, &dofs, ref_cell, &fe_trafo, ref_elem, modal_basis, conform);
  }
  else if (fe_type == FE_TYPE_BDM)
  {  
    this->lagrange_only_ = false;
    assert (param.size() == 1);
    assert (param[0] > 0);
    create_BDM_element(topo_cell_type, param[0], &ansatz, &dofs, ref_cell, &fe_trafo, ref_elem, modal_basis, conform);
  }
  else if (fe_type == FE_TYPE_RT)
  { 
    this->lagrange_only_ = false;
    assert (param.size() == 1);
    assert (param[0] >= 0);
    create_RT_element(topo_cell_type, param[0], &ansatz, &dofs, ref_cell, &fe_trafo, ref_elem, modal_basis, conform);
  }
  else
  {
    std::cout << "Unexpected fe type " << std::endl;
    quit_program();
  }
  // initialize reference element
  assert (ansatz != nullptr);
  assert (dofs != nullptr);
  assert (fe_trafo != nullptr);
  assert (ref_cell);
  
  assert (ansatz->type() != ANSATZ_NOT_SET);
  assert (dofs->type() != DOF_CONTAINER_NOT_SET);
  
  ref_elem->init(ansatz, dofs, fe_trafo, modal_basis, fe_type);
  ref_elem->set_instance_id(this->ref_elements_.size());
  this->ref_elements_.push_back(ref_elem);
  this->ref_cells_.push_back(ref_cell);
  this->dof_containers_.push_back(dofs);
  this->fe_trafos_.push_back(fe_trafo);
  this->ansatz_spaces_.push_back(ansatz);
  this->max_fe_conform_.push_back(conform);
  
  this->added_fe_types_.push_back(fe_type);
  this->added_cell_types_.push_back(topo_cell_type);
  this->added_params_.push_back(param);
  
  return this->ref_elements_.size()-1;
}

template <class DataType, int DIM>
CellTrafoPtr<DataType, DIM> FEInstance<DataType, DIM>::create_cell_trafo (size_t fe_id, std::vector<DataType> coord_vtx) const
//CellTrafoPtr<DataType, DIM> FEInstance<DataType, DIM>::create_cell_trafo (size_t fe_id, int align_number) const
{
  assert (fe_id < this->ref_elements_.size());
  ConstRefCellPtr<DataType, DIM> ref_cell = this->ref_cells_[fe_id];
  
  CellTrafoPtr<DataType, DIM> cell_trafo = nullptr;

  switch(ref_cell->type())
  {
    case REF_CELL_LINE_STD:
      cell_trafo = CellTrafoPtr<DataType, DIM> (new LinearLineTransformation<DataType, DIM>(ref_cell));
      break;
    case REF_CELL_TRI_STD:
      cell_trafo = CellTrafoPtr<DataType, DIM> (new LinearTriangleTransformation<DataType, DIM>(ref_cell));
      break;
    case REF_CELL_TET_STD:
      cell_trafo = CellTrafoPtr<DataType, DIM> (new LinearTetrahedronTransformation<DataType, DIM>(ref_cell));
      break;
    case REF_CELL_PYR_STD:
      cell_trafo = CellTrafoPtr<DataType, DIM> (new LinearPyramidTransformation<DataType, DIM>(ref_cell));
      break;
    case REF_CELL_QUAD_STD:
#ifdef FORCE_NONLINEAR_TRAFO
      if (mesh::is_parallelogram(coord_vtx) && false)
#else
      if (mesh::is_parallelogram(coord_vtx))
#endif
      {
        cell_trafo = CellTrafoPtr<DataType, DIM> (new AlignedQuadTransformation<DataType, DIM>(ref_cell));
      }
      else
      {
        cell_trafo = CellTrafoPtr<DataType, DIM> (new BiLinearQuadTransformation<DataType, DIM>(ref_cell));
        this->nb_nonlin_trafos_++;
      }
      break;
    case REF_CELL_HEX_STD:
#ifdef FORCE_NONLINEAR_TRAFO
      if (mesh::is_parallelepiped(coord_vtx) && false)
#else
      if (mesh::is_parallelepiped(coord_vtx))
#endif
      {
        cell_trafo = CellTrafoPtr<DataType, DIM> (new AlignedHexahedronTransformation<DataType, DIM>(ref_cell));
      }
      else
      {
        cell_trafo = CellTrafoPtr<DataType, DIM> (new TriLinearHexahedronTransformation <DataType, DIM>(ref_cell));
        this->nb_nonlin_trafos_++;
      }
      break;
    default:
      std::cout << "Unexpected ref cell type " << std::endl;
      quit_program();
  }
  cell_trafo->reinit(coord_vtx);
  return cell_trafo;
}

template class FEInstance< float, 3 >;
template class FEInstance< float, 2 >;
template class FEInstance< float, 1 >;

template class FEInstance< double, 3 >;
template class FEInstance< double, 2 >;
template class FEInstance< double, 1 >;

} // namespace doffem
} // namespace hiflow
