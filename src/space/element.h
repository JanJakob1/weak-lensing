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

#ifndef HIFLOW_ELEMENT_H
#define HIFLOW_ELEMENT_H

#include "common/vector_algebra.h"
#include "common/log.h"
#include "dof/dof_fem_types.h"
#include "mesh/types.h"
#include "mesh/entity.h"

/// @author Staffan Ronnas, Philipp Gerstner

namespace hiflow {

namespace mesh {
class Mesh;
//class Entity;
}

namespace doffem {
template <class DataType, int DIM> class RefElement;
template <class DataType, int DIM> class FEManager;
template <class DataType, int DIM> class CellTransformation;
template <class DataType, int DIM> class RefCell;
}

template < class T > class Quadrature;
template <class DataType, int DIM> class VectorSpace;

/// \brief Class representing one (physical) element of a finite element space. 
/// In case of vector valued spaces consisting of tensor products of different finite elements, 
/// this class acts as a vector valued function space.
///
/// Important: it is assumed that all individuals RefElements live on the same reference cell
///
/// \details This class provides a local view of one element of
/// the finite element space, and brings together information from
/// the mesh, dof and fem data structures.

template < class DataType, int DIM > 
class Element 
{
public:
  typedef Vec<DIM, DataType> Coord;
  typedef boost::function< size_t (size_t, size_t) > IndexFunction;

  /// \brief Construct an element on a given cell in a space.
  ///
  /// \param[in] space         finite element space to which element belongs
  /// \param[in] cell_index    index of cell for which element should be created

  Element(const VectorSpace< DataType, DIM > &space, int cell_index);

  ~Element()
  {}

  /////////////////////////// Simple informations ////////////////////////////////////////
  
  /// \return the number of variables in the finite element space
  size_t nb_var() const;

  /// \return the number of finite elements in the finite element space.
  /// This value might differ from nb_var() if vector valued elements are used that are 
  /// not obtained by tensor product of single elements.
  size_t nb_fe() const;

  size_t var_2_fe (size_t var) const;

  size_t var_2_comp (size_t var) const;

  std::vector<size_t> fe_2_var(size_t fe_ind) const;
  
  /// \brief Access the number of dofs for a given variable associated with the
  /// element. 
  /// \param[in]  var   the number of the variable 
  /// \return the number of dofs associated with variable var
  size_t nb_dof(size_t fe_ind) const;

  size_t nb_comp(size_t fe_ind) const;
  
  size_t dof_offset (size_t fe_ind) const;
  
  /// \return the cell index of the element.
  int cell_index() const ;
  
  /// \return true if and only if the element belongs to the local subdomain.
  bool is_local() const;
  
  /// \return true if and only if the element is adjacent to the boundary of the
  /// mesh.
  bool is_boundary() const;

  /// \return the local facet numbers of the element which lie on the boundary.
  std::vector< int > boundary_facet_numbers() const;
  
  /////////////////////////// Get pointers to complex objects ////////////////////////////
  /// \return a reference to the finite element space to which the element belongs.
  inline const VectorSpace< DataType, DIM > &get_space() const 
  { 
    return *this->space_; 
  }

  /// \return a reference to the cell entity associated with the element.
  inline const mesh::Entity &get_cell() const 
  { 
    return this->cell_;
  }

  /// \return the cell transformation associated with the element
  std::shared_ptr<doffem::CellTransformation< DataType, DIM > > get_cell_transformation(/*size_t fe_ind=0*/) const;

  /// \brief Accesss the finite element ansatz for a variable on the cell.
  /// \param[in] var    the number of the variable
  /// \return a pointer to the finite element ansatz for variable var on the
  /// cell.
  const doffem::RefElement< DataType, DIM > *get_fe(size_t fe_ind) const;

  const doffem::RefElement< DataType, DIM > *get_fe_for_var(size_t var) const;

  doffem::ConstRefCellPtr<DataType, DIM> ref_cell() const;

  /////////////////////////// Dof index handling /////////////////////////////////////////
  /// \brief Access the dof indices for a given variables associated with the
  /// element. \param[in]  var       number of the variable \param[out] indices
  /// vector of dof indices for variable @p var associated with the element.
  void get_dof_indices(size_t fe_ind, std::vector< int > &indices) const;

  /// \brief Access the dof indices associated with the element.
  /// \param[out] indices   vector of dof indices associated with the element.
  void get_dof_indices(std::vector< int > &indices) const;

  /// \brief Access the dof indices on the boundary for a given variable
  /// associated with the element. \param[in]  var       number of the variable
  /// \param[out] indices   vector of dof indices for variable @p var associated
  /// with the element.
  void get_dof_indices_on_subentity(size_t fe_ind, int tdim, int sindex,
                                    std::vector< int > &indices) const;

  /// \brief Access the dof indices on the boundary associated with the element.
  /// \param[out] indices   vector of dof indices associated with the element.
  void get_dof_indices_on_subentity(int tdim, int sindex,
                                    std::vector< int > &indices) const;

  /////////////////////////// Evaluation of basis functions //////////////////////////////
  /// indexing of return array of routines N(pt, weight), grad_N(ot, weight) and hessians_N(pt, weight)
  /// @param[in] var considered (physical) variable
  /// @param[in] i index of basis function in function space associated to given variable
  /// Note: 0 <= i < nb_dof( var_2_fe(var) ) 
  /// \return index in array
  size_t iv2ind (size_t i, size_t var) const;

  /// evaluate all components of all (mapped) basis functions associated to specified finite element of index fe_ind
  /// -> return.size = sum_{refelements fe} fe.nb_comp * fe.dim
  /// indexing the return array is possible with routine iv2ind()
  /// important: for performance reasons, the user has tp provide the coordinate on the 
  /// corresponding reference cell
  
  void N(const Coord &ref_pt,   
         std::vector< DataType > &vals) const;
         
  void grad_N(const Coord &ref_pt, 
              std::vector< Vec<DIM,DataType> > &gradients) const;
              
  void hessian_N(const Coord &ref_pt, 
                 std::vector< Mat<DIM, DIM, DataType> > &hessians) const;

  /// evaluate all components of all (mapped) basis functions associated to specified finite element of index fe_ind
  /// -> return.size = fe(fe_ind).nb_comp * fe(fe_ind).dim
  /// indexing the return array is possible with routine get_fe(fe_ind)->iv2ind()
  /// important: for performance reasons, the user has tp provide the coordinate on the 
  /// corresponding reference cell

  void N_fe(const Coord &ref_pt, 
            size_t fe_ind, 
            std::vector< DataType > &vals) const;
  
  void grad_N_fe(const Coord &ref_pt, 
                 size_t fe_ind, 
                 std::vector< Vec<DIM,DataType> > &gradients) const;
  
  void N_and_grad_N_fe (const Coord &ref_pt, 
                        size_t fe_ind, 
                        std::vector< DataType > &vals,
                        std::vector< Vec<DIM, DataType> > &gradients) const;
                                              
  void hessian_N_fe(const Coord &ref_pt, 
                    size_t fe_ind, 
                    std::vector< Mat<DIM, DIM, DataType> > &hessians) const;

  /// evaluate all (mapped) basis functions associated to specified variable
  /// -> return.size = get_fe(var_2_fe(var)).dim
  /// important: for performance reasons, the user has tp provide the coordinate on the 
  /// corresponding reference cell

  void N_var(const Coord &ref_pt, 
             size_t var, 
             std::vector< DataType > &vals) const;
             
  void grad_N_var(const Coord &ref_pt, 
                  size_t var, 
                  std::vector< Vec<DIM,DataType> > &gradients) const;
                  
  void hessian_N_var(const Coord &ref_pt, 
                     size_t var, 
                     std::vector< Mat<DIM, DIM, DataType> > &hessians) const;
private:

  std::vector<size_t> weight_offsets_;
  std::vector<size_t> dim_offsets_;
  size_t nb_comp_;
  size_t weight_size_;
  size_t dim_;
  
  size_t active_fe_ind_;
  std::vector<DataType> active_dof_values_;
  bool mapping_eval_;
  
  const VectorSpace< DataType, DIM > *space_;
  mesh::Entity cell_;
};



} // namespace hiflow

#endif
