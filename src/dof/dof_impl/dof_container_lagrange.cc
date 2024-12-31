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

#include "common/data_tools.h"
#include "dof/dof_impl/dof_container_lagrange.h"
#include "dof/dof_impl/dof_functional_point.h"
#include "mesh/geometric_tools.h"
#include "mesh/entity.h"

namespace hiflow {
namespace doffem {

template<class DataType, int DIM>
DofContainerLagrange<DataType, DIM>::~DofContainerLagrange()
{
  this->clear();
}

template<class DataType, int DIM>
void DofContainerLagrange<DataType, DIM>::evaluate (FunctionSpace<DataType, DIM> const * space, 
                                                    const std::vector< DofID > & dof_ids, 
                                                    std::vector< std::vector<DataType> >& dof_values ) const
{
  assert (this->initialized_);
  assert (space != nullptr);

  dof_values.resize(dof_ids.size());
  
  // loop over dofs
  for (size_t l = 0; l < dof_ids.size(); ++l)
  {
    DofID dof_id = dof_ids[l];
    assert (dof_id < this->dofs_.size());
    assert (this->dofs_[dof_id] != nullptr);
    
    set_to_zero(space->dim(), dof_values[l]);
        
    dynamic_cast< DofPointEvaluation<DataType, DIM> const * >(this->dofs_[dof_id])->evaluate (space, dof_values[l]);
  }
}

template < class DataType, int DIM > 
void DofContainerLagrange<DataType, DIM>::evaluate (RefCellFunction<DataType, DIM> const * func, 
                                                    const std::vector< DofID > & dof_ids, 
                                                    std::vector< std::vector<DataType> >& dof_values ) const 
{
  assert (this->initialized_);
  assert (func != nullptr);
  
  const size_t nb_func = func->nb_func();
  const size_t nb_comp = func->nb_comp();
   
  dof_values.resize(dof_ids.size());
  
  // loop over all specified dof functionals
  for (size_t l = 0; l < dof_ids.size(); ++l)
  {
    const DofID dof_id = dof_ids[l];
    assert (dof_id >= 0);
    assert (dof_id < this->dim_);
    assert (this->dofs_[dof_id] != nullptr);
    
    set_to_zero(nb_func, dof_values[l]);
   
    dynamic_cast< DofPointEvaluation<DataType, DIM> const * >(this->dofs_[dof_id])->evaluate (func, 0, dof_values[l]); 
  }
}

template<class DataType, int DIM>
std::vector<Vec<DIM, DataType> > DofContainerLagrange<DataType, DIM>::get_dof_coords() const
{
  std::vector< Vec<DIM, DataType> > coords;
  for (size_t l=0; l<this->dofs_.size(); ++l)
  {
    DofPointEvaluation<DataType, DIM> const * cur_dof = 
      dynamic_cast< DofPointEvaluation<DataType, DIM> const * >(this->dofs_[l]);

    assert (cur_dof != 0);

    coords.push_back(cur_dof->get_point());
  }
  return coords;
}
template<class DataType, int DIM>
std::vector< Vec<DIM, DataType> > DofContainerLagrange<DataType, DIM>::get_dof_coords_on_subentity (size_t tdim, int sindex) const
{
  std::vector<DofID> dof_ids_on_subentity = this->get_dof_on_subentity(tdim, sindex);
  
  std::vector< Vec<DIM, DataType> > coords;
  for (size_t l=0; l<dof_ids_on_subentity.size(); ++l)
  {
    DofPointEvaluation<DataType, DIM> const * cur_dof = 
      dynamic_cast< DofPointEvaluation<DataType, DIM> const * >(this->dofs_[dof_ids_on_subentity[l]]);

    assert (cur_dof != 0);

    coords.push_back(cur_dof->get_point());
  }
  return coords;
  
  
  
}

template<class DataType, int DIM>
size_t DofContainerLagrange<DataType, DIM>::get_degree() const
{
  return degree_;
}

template<class DataType, int DIM>
void DofContainerLagrange<DataType, DIM>::init (size_t degree, size_t nb_comp)
{
  if (this->initialized_ && this->degree_ == degree && this->nb_comp_ == nb_comp)
  {
    return;
  }
  
  this->clear();
  
  std::vector<Coord> dof_coords;
  this->degree_ = degree;
  this->nb_comp_ = nb_comp;
  
  RefCellType ref_cell_type = this->ref_cell_->type();
  // initialize reference cell and dof points
  switch (ref_cell_type)
  {
    case REF_CELL_LINE_STD:
      this->init_dof_coord_line(degree, dof_coords);
      break;
    case REF_CELL_TRI_STD:
      this->init_dof_coord_tri(degree, dof_coords);
      break;
    case REF_CELL_QUAD_STD:
      this->init_dof_coord_quad(degree, dof_coords);
      break;
    case REF_CELL_TET_STD:
      this->init_dof_coord_tet(degree, dof_coords);
      break;
    case REF_CELL_HEX_STD:
      this->init_dof_coord_hex(degree, dof_coords);
      break;
    case REF_CELL_PYR_STD:
      this->init_dof_coord_pyr(degree, dof_coords);
      break;
    default:
      std::cerr << "Unexpected reference cell type " << std::endl;
      quit_program();
  }
    
  // initialize dof functionals
  const size_t nb_dof_coord = dof_coords.size();
  for (size_t c=0; c<this->nb_comp_; ++c)
  { 
    for (size_t i=0; i<nb_dof_coord; ++i)
    {
      DofPointEvaluation<DataType, DIM> * dof = new DofPointEvaluation<DataType, DIM>();
      dof->init (dof_coords[i], c, this->ref_cell_);
      this->push_back( dof );
    }
  }
  DofContainer<DataType, DIM>::init();
}

template<class DataType, int DIM>
void DofContainerLagrange<DataType, DIM>::init_dof_coord_line (size_t deg, std::vector<Coord> &dof_coords) const
{
  assert (DIM == 1);
  dof_coords.clear();

  if (deg == 0) 
  {
    // Coordinates of the middle-point of line
    Coord coord;

    // Centroid
    coord[0] = 0.5;

    dof_coords.push_back(coord);
  } 
  else 
  {
    assert(deg > 0);

    // Lexicographical ordering
    const DataType offset = (1.0 / deg);
    const size_t nb_dof_on_cell = (deg + 1);
    dof_coords.resize(nb_dof_on_cell);

    const size_t nb_dof_on_line = deg + 1;

    // Filling coord vector by lexicographical strategy
    for (size_t j = 0; j < nb_dof_on_line; ++j) 
    {
      Coord coord;
      coord[0] = j * offset;

      dof_coords[j] = coord;
    }
  }
}

template<class DataType, int DIM>
void DofContainerLagrange<DataType, DIM>::init_dof_coord_tri (size_t deg, std::vector<Coord> &dof_coords) const
{
  assert (DIM == 2);
  dof_coords.clear();
  
  // Lexicographical ordering
  if (deg == 0) 
  {
    // Coordinates of the middle-point of triangle
    Coord coord;

    // Centroid
    coord[0] = 1.0 / 3.0;
    coord[1] = 1.0 / 3.0;

    dof_coords.push_back(coord);
  } 
  else 
  {
    assert(deg > 0);
    DataType offset = (1.0 / deg);

    size_t nb_dof_on_cell = (((deg + 1) * (deg + 2)) / 2);
    dof_coords.resize(nb_dof_on_cell);

    // Filling coord vector for full triangle by lexicographical strategy
    // and filtering the line coordinates by given mesh ordering strategy
    // which was computed in first step

    const size_t nb_dof_on_line = deg + 1;

    for (size_t j = 0; j < nb_dof_on_line; ++j) // y axis
    { 
      const DataType j_offset = j * offset;
      
      size_t j_ind_offset = 0;
      for (size_t n = 0; n < j; ++n)
      {
        j_ind_offset += nb_dof_on_line - n;
      }
      for (size_t i = 0; i < nb_dof_on_line - j; ++i) // x axis
      {
        Coord coord;

        coord[0] = i * offset;
        coord[1] = j_offset;

        dof_coords[i + j_ind_offset] = coord;
//      std::cout << coord[0] << " , " << coord[1] << std::endl; 
      }
    }
  }
}

template<class DataType, int DIM>
void DofContainerLagrange<DataType, DIM>::init_dof_coord_quad (size_t deg, std::vector<Coord> &dof_coords) const
{
  assert (DIM == 2);
  dof_coords.clear();

  if (deg == 0) 
  {
    // Coordinates of the middle-point of quad
    Coord coord;

    // Centroid
    coord[0] = 0.5;
    coord[1] = 0.5;

    dof_coords.push_back(coord);
  } 
  else 
  {
    assert(deg > 0);

    // Lexicographical ordering

    const DataType offset = (1.0 / deg);
    const size_t nb_dof_on_cell = (deg + 1) * (deg + 1);
    dof_coords.resize(nb_dof_on_cell);

    const size_t nb_dof_on_line = deg + 1;


    // Filling coord vector by lexicographical strategy
    for (size_t j = 0; j < nb_dof_on_line; ++j) 
    {
      const DataType j_offset = j * offset;
      for (size_t i = 0; i < nb_dof_on_line; ++i) 
      {
        Coord coord;

        coord[0] = i * offset;
        coord[1] = j_offset;

        dof_coords[i + j * nb_dof_on_line] = coord;
      }
    }
  }
}

template<class DataType, int DIM>
void DofContainerLagrange<DataType, DIM>::init_dof_coord_tet (size_t deg, std::vector<Coord> &dof_coords) const
{
  assert (DIM == 3);
  dof_coords.clear();

  // Lexicographical ordering

  if (deg == 0) 
  {
    // Coordinates of the middle-point of tet
    Coord coord;

    // Centroid
    coord[0] = 0.25;
    coord[1] = 0.25;
    coord[2] = 0.25;

    dof_coords.push_back(coord);
  } 
  else 
  {
    assert(deg > 0);

    const DataType offset = (1.0 / deg);

    const size_t nb_dof_on_cell = (deg + 1) * (deg + 2) * (deg + 3) / 6;
    dof_coords.resize(nb_dof_on_cell);

    // Filling coord vector for full tet by lexicographical strategy
    // and filtering the line coordinates by given mesh ordering strategy
    // which was computed in first step

    const size_t nb_dof_on_line = deg + 1;  

    for (size_t k = 0; k < nb_dof_on_line; ++k) 
    { // z axis
      const DataType k_offset = k * offset;
      
      // First: index offset z axis
      size_t k_ind_offset = 0;
      for (size_t m = 0; m < k; ++m) 
      {
        const size_t help = nb_dof_on_line - m;
        for (size_t dof = 0; dof < nb_dof_on_line - m; ++dof) 
        {
          k_ind_offset += help - dof;
        }
      }

      for (size_t j = 0; j < nb_dof_on_line - k; ++j) 
      { // y axis
        const DataType j_offset = j * offset;

         // Second: increasing index offset by y axis on current z axis
        size_t ind_offset = k_ind_offset;
        for (size_t n = 0; n < j; ++n) 
        {
          ind_offset += nb_dof_on_line - n - k;
        }

        for (size_t i = 0; i < nb_dof_on_line - k - j; ++i) // x axis
        {
          Coord coord;
          coord[0] = i * offset;
          coord[1] = j_offset;
          coord[2] = k_offset;

          dof_coords[i + ind_offset] = coord;
        }
      }
    }
  }
}

template<class DataType, int DIM>
void DofContainerLagrange<DataType, DIM>::init_dof_coord_hex (size_t deg, std::vector<Coord> &dof_coords) const
{
  assert (DIM == 3);
  dof_coords.clear();
  
  // set coordinates of all DoF points in lexicographical ordering
  if (deg == 0) 
  {
    // Coordinates of middle-point of hex
    Coord coord;

    coord[0] = 0.5;
    coord[1] = 0.5;
    coord[2] = 0.5;

    dof_coords.push_back(coord);
  } 
  else 
  {
    assert(deg > 0);

    const DataType offset = (1.0 / deg);
    const size_t fd_1 = deg + 1;

    const size_t nb_dof_on_cell = fd_1 * fd_1 * fd_1;
    const size_t nb_dof_on_line = fd_1;

    dof_coords.resize(nb_dof_on_cell);

    /*     Triquadratic Hexaedron (Q2-Element)
           ===================================

                24----------25-----------26
                /|          /|          /|
               / |         / |         / |
             15----------16-----------17 |
             /|  |       /|  |       /|  |
            / | 21------/---22------/-|--23
           6-----------7-----------8  | /|
           |  |/ |     |  |/ |     |  |/ |
           | 12--|-----|-13--|-----|--14 |
           | /|  |     | /|  |     | /|  |
           |/ | 18-----|----19-----|/-|--20
           3-----------4-----------5  | /
           |  |/       |  |/       |  |/
           |  9--------|-10--------|--11
           | /         | /         | /
           |/          |/          |/
           0-----------1-----------2          */

    // Filling coord vector for full hex by lexicographical strategy

    for (size_t k = 0; k < nb_dof_on_line; ++k) 
    {
      const DataType k_offset = k * offset;
      for (size_t j = 0; j < nb_dof_on_line; ++j) 
      {
        const DataType j_offset = j * offset;
        for (size_t i = 0; i < nb_dof_on_line; ++i) 
        {
          Coord coord;
          coord[0] = i * offset;
          coord[1] = j_offset;
          coord[2] = k_offset;

          dof_coords[(i + (j + k * nb_dof_on_line) * nb_dof_on_line)] = coord;
        }
      }
    }
  }
}

template<class DataType, int DIM>
void DofContainerLagrange<DataType, DIM>::init_dof_coord_pyr (size_t deg, std::vector<Coord> &dof_coords) const
{
  assert (deg <= 2);
  assert (DIM == 3);

  dof_coords.clear();

  // Lexicographical ordering
  if (deg == 0) 
  {
    Coord coord;

    // Centroid
    coord[0] = 0.5;
    coord[1] = 0.5;
    coord[2] = 0.25;

    dof_coords.push_back(coord);
  } 
  else 
  {
    assert(deg > 0);

    const DataType offset = (1.0 / deg);
    size_t nb_dof_on_cell;

    if (deg == 1) 
    {
      nb_dof_on_cell = 5;
    } 
    else if (deg == 2)
    {
      nb_dof_on_cell = 14;
    } 
    else 
    {
      std::cerr << "Unexpected finite element degree " << deg << std::endl;
      quit_program();
    }

    dof_coords.resize(nb_dof_on_cell);

    const size_t nb_dof_on_line = deg + 1;

    for (size_t k = 0; k < nb_dof_on_line; ++k) // z axis
    {       
      size_t k_ind_offset = 0;
      for (size_t m = 0; m < k; ++m) 
      {
        const size_t help = nb_dof_on_line - m;
        k_ind_offset += help * help;
      }

      for (size_t j = 0; j < nb_dof_on_line - k; ++j) // y axis 
      { 
        size_t j_ind_offset = k_ind_offset;
        for (size_t n = 0; n < j; ++n) 
        {
          j_ind_offset += nb_dof_on_line - k;
        }

        for (size_t i = 0; i < nb_dof_on_line - k; ++i) // x axis
        {
          Coord coord;

          coord[0] = k * offset * 0.5 + i * offset;
          coord[1] = k * offset * 0.5 + j * offset;
          coord[2] = k * offset;

          dof_coords[i + j_ind_offset] = coord;
        }
      }
    }
  }
}

template class DofContainerLagrange< float, 3 >;
template class DofContainerLagrange< float, 2 >;
template class DofContainerLagrange< float, 1 >;

template class DofContainerLagrange< double, 3 >;
template class DofContainerLagrange< double, 2 >;
template class DofContainerLagrange< double, 1 >;

} // namespace doffem
} // namespace hiflow
