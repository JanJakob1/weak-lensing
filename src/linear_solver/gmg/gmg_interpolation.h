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

/// \author Aksel Alpay, Martin Wlotzka

#ifndef GMG_INTERPOLATION_H
#define GMG_INTERPOLATION_H

#include "dof_identification.h"

namespace hiflow {
namespace la {
namespace gmg {

typedef std::vector< int > ijk_dof_id;

/// ijk iterators traverse the dofs of a cell type in the same order
/// as the local dof ids are numbered. They store the current local dof id
/// and the corresponding ijk id, and can thus be used to establish
/// a relationship between ijk and local dof ids

class IjkIterator {
public:
  /// Construct ijk iterator
  /// @param dimension The topological dimension of the cell
  /// @param fe_degree The degree of the finite element ansatz

  IjkIterator(unsigned dimension, unsigned fe_degree)
      : current_ijk_id_(dimension, 0), current_local_id_(0),
        fe_degree_(fe_degree), dimension_(dimension), is_end_(false) {
    assert(dimension > 0);
    // currently more than three dimensions are not supported
    assert(dimension <= 3);

    assert(fe_degree > 0);
  }

  virtual ~IjkIterator() {}

  /// move to the next dof
  virtual void next() = 0;

  /// @return whether the iterator can be incremented further. true will be
  /// returned one step after the last dof has been reached.

  bool is_at_end() const { return is_end_; }

  /// resets the iterator to the first dof of the cell

  void reset() {
    current_local_id_ = 0;
    current_ijk_id_ = ijk_dof_id(current_ijk_id_.size(), 0);
    is_end_ = false;
  }

  /// @return the ijk dof id at the current position

  const ijk_dof_id &get_current_ijk_id() const { return current_ijk_id_; }

  /// @return the (cell-) local dof id at the current position

  int get_current_local_id() const { return current_local_id_; }

  /// @return the dimension of the cell

  unsigned get_dimension() const { return dimension_; }

  /// @return the degree of the fe ansatz

  unsigned get_fe_degree() const { return fe_degree_; }

  /// move to the next dof
  /// @return A reference to this ijk iterator object

  IjkIterator &operator++() {
    next();
    return *this;
  }

  /// move to the next dof.
  /// To avoid unboxing issues with derived classes, the postincrement operator
  /// behaves the same way the preincrement operator does.
  /// @return A reference to this ijk iterator object

  IjkIterator &operator++(int) {
    ++(*this);
    return *this;
  }

  /// advance the iterator several steps
  /// @param n The number of steps to advance
  /// @return A reference to this ijk iterator object

  IjkIterator &operator+=(int n) {
    for (; n > 0; --n)
      next();

    return *this;
  }

protected:
  bool is_end_;

  ijk_dof_id current_ijk_id_;
  int current_local_id_;
  const unsigned fe_degree_;
  const unsigned dimension_;
};

/// Provides an ijk iterator implentation for quadrilaterals and hexahedrons

class QuadBasedIjkIterator : public IjkIterator {
public:
  /// @param dimension The topological dimension of the cell. Only two
  /// and three dimensions are supported.
  /// @param fe_degree The fe ansatz

  QuadBasedIjkIterator(unsigned dimension, unsigned fe_degree)
      : IjkIterator(dimension, fe_degree), is_three_dimensional_(false) {
    assert(dimension == 2 || dimension == 3);

    if (dimension == 3)
      is_three_dimensional_ = true;
  }

  virtual ~QuadBasedIjkIterator() {}

  /// advance the iterator to the next dof

  virtual void next() {
    if (is_end_)
      return;

    if (current_ijk_id_[0] < fe_degree_)
      ++(current_ijk_id_[0]);
    else {
      if (current_ijk_id_[1] < fe_degree_) {
        ++(current_ijk_id_[1]);
        current_ijk_id_[0] = 0;
      } else {
        if (is_three_dimensional_) {
          if (current_ijk_id_[2] < fe_degree_) {
            ++(current_ijk_id_[2]);
            current_ijk_id_[0] = 0;
            current_ijk_id_[1] = 0;
          } else {
            is_end_ = true;
            return;
          }
        } else {
          // We have reached the final dof, and cannot go further.
          is_end_ = true;
          return;
        }
      }
  }

  ++current_local_id_;
}

private:
  bool is_three_dimensional_;
};

/// Implements an ijk iterator for quadrilaterals

class QuadrilateralIjkIterator : public QuadBasedIjkIterator {
public:
  /// @param fe_degree The degree of the fe ansatz

  explicit QuadrilateralIjkIterator(unsigned fe_degree)
      : QuadBasedIjkIterator(2, fe_degree) {}

  virtual ~QuadrilateralIjkIterator() {}
};

/// Implements an ijk iterator for hexahedrons

class HexahedronIjkIterator : public QuadBasedIjkIterator {
public:
  /// @param fe_degree The degree of the fe ansatz

  explicit HexahedronIjkIterator(unsigned fe_degree)
      : QuadBasedIjkIterator(3, fe_degree) {}

  virtual ~HexahedronIjkIterator() {}
};

/// Implements an ijk iterator for triangles

class TriangleIjkIterator : public IjkIterator {
public:
  /// @param fe_degree The degree of the fe ansatz

  explicit TriangleIjkIterator(unsigned fe_degree)
      : IjkIterator(2, fe_degree) {}

  virtual ~TriangleIjkIterator() {}

  /// advance to the next dof

  virtual void next() {
    if (is_end_)
      return;

    if (current_ijk_id_[0] < last_index_in_row(current_ijk_id_[1]))
      ++(current_ijk_id_[0]);
    else {
      if (current_ijk_id_[1] < fe_degree_) {
        ++(current_ijk_id_[1]);
        current_ijk_id_[0] = 0;
      } else {
        // We have reached the final dof, and cannot go further.
        is_end_ = true;
        return;
      }
    }

    ++current_local_id_;
  }

private:
  /// @return In ijk notation: the last valid i index in a row j
  /// @param row_index the j index of the row

  inline unsigned last_index_in_row(int row_index) const {
    return fe_degree_ - row_index;
  }
};

/// Implements an ijk iterator for tetrahedrons

class TetrahedronIjkIterator : public IjkIterator {
public:
  /// @param fe_degree The degree of the fe ansatz

  explicit TetrahedronIjkIterator(unsigned fe_degree)
      : IjkIterator(3, fe_degree) {}

  virtual ~TetrahedronIjkIterator() {}

  /// advance to the next dof

  virtual void next() {
    if (is_end_)
      return;

    if (current_ijk_id_[0] <
        last_index_in_row(current_ijk_id_[1], current_ijk_id_[2]))
      ++(current_ijk_id_[0]);
    else {
      if (current_ijk_id_[1] < last_index_in_layer(current_ijk_id_[2])) {
        ++(current_ijk_id_[1]);
        current_ijk_id_[0] = 0;
      } else {
        if (current_ijk_id_[2] < fe_degree_) {
          ++(current_ijk_id_[2]);
          current_ijk_id_[0] = 0;
          current_ijk_id_[1] = 0;
        } else {
          // We have reached the final dof, and cannot go further.
          is_end_ = true;
          return;
        }
      }
    }

    ++current_local_id_;
  }

private:
  /// @return In ijk notation: the last valid j index in a layer k
  /// @param layer_index the k index of the layer

  inline unsigned last_index_in_layer(int layer_index) const {
    return fe_degree_ - layer_index;
  }

  /// @return In ijk notation: the last valid i index in a row j
  /// @param row_index the j index of the row
  /// @param layer_index the k index of the row

  inline unsigned last_index_in_row(int row_index, int layer_index) const {
    return last_index_in_layer(layer_index) - row_index;
  }
};

/// For a given cell and variable, returns an ijk iterator of corresponding type

template < class LAD, int DIM > class IjkIteratorFactory {
public:
  TYPE_FROM_CLASS(LAD, DataType);
  /// @param element The element for which to construct ijk iterators

  explicit IjkIteratorFactory(const Element< DataType, DIM > *element)
      : element_(element) {}

  /// @return an ijk iterator object that matches the cell and variable.
  /// memory management is the user's responsibility
  /// @param var The variable for which to construct an ijk iterator

  IjkIterator *get(int var) const {
    CellType::Tag t = element_->get_cell().cell_type().tag();


  
    DofContainerLagrange <DataType, DIM > const * dof_element_ = 
        dynamic_cast < DofContainerLagrange <DataType, DIM > const * >(element_->get_fe(var)->dof_container());

    int degree =  dof_element_->get_degree();


    switch (t) {
    case CellType::QUADRILATERAL:
      return new QuadrilateralIjkIterator(degree);
    case CellType::TRIANGLE:
      return new TriangleIjkIterator(degree);
    case CellType::HEXAHEDRON:
      return new HexahedronIjkIterator(degree);
    case CellType::TETRAHEDRON:
      return new TetrahedronIjkIterator(degree);
    default:
      throw std::invalid_argument("invalid cell type in IjkIteratorFactory");
    }
  }

private:
  const Element< DataType, DIM > *element_;
};


/// Provides methods to translate global dof ids (relative to the whole mesh),
/// local dof ids (relative to the cell) and ijk dof ids (relative to the cell
/// and variable) into each other.

template < class LAD, int DIM > 
class DofIdConverter {
  TYPE_FROM_CLASS(LAD, DataType);

  UMAP_NAMESPACE_GMG::unordered_map< int, int > g2lmap_;

  /// Translates a local dof id to its ijk representation.
  /// local_to_ijk_map_[var][i] = ijk representation of local DoF with id i
  std::vector< std::vector< ijk_dof_id > > l2ijk_map_;

  typedef util::multi_array< int > IjkToLocalMapType;

  std::vector< IjkToLocalMapType > ijk2l_map_;

public:
  /// @param element The element for which dofs shall be translated

  explicit DofIdConverter(const Element< DataType, DIM > *element) {
    assert(element);

    tag_ = element->get_cell().cell_type().tag();

    // build global -> local map

    int total_num_dofs = 0;
    for (int var = 0; var < element_->nb_fe(); ++var)			//TODO: nb_fe correct?
      total_num_dofs += element_->nb_dof(var);

    g2lmap_ = UMAP_NAMESPACE_GMG::unordered_map< int, int >(1.2 * total_num_dofs);

    init_g2l_map();

    // next, build ijk translation data structures

    l2ijk_map_.resize(element->nb_fe(), std::vector< ijk_dof_id >());		//TODO: nb_fe correct?
    ijk2l_map_.resize(element->nb_fe());

    for (int var = 0; var < element->nb_fe(); ++var) {
      IjkIteratorFactory< LAD, DIM > iter_fact(element);
      boost::shared_ptr< IjkIterator > ijk_iter(iter_fact.get(var));

      int num_dofs = element->nb_dof(var);

      l2ijk_map_[var] = std::vector< ijk_dof_id >(
          num_dofs, ijk_dof_id(ijk_iter->get_dimension()));

      unsigned fe_degree = ijk_iter->get_fe_degree();

      unsigned k_size = 1;
      if (ijk_iter->get_dimension() == 3) {
        util::multi_array< int >::size_type array_sizes[] = {
            fe_degree + 1, fe_degree + 1, fe_degree + 1};

        ijk2l_map_[var] = IjkToLocalMapType(array_sizes);
      } else {
        util::multi_array< int >::size_type array_sizes[] = {fe_degree + 1,
                                                             fe_degree + 1};

        ijk2l_map_[var] = IjkToLocalMapType(array_sizes);
      }

      // Unused values of the ijk->local map are supposed to be -1,
      // thus we fill the map with -1 values.
      // This is important because the ijk2l_map_ will, as a multidimensional
      // array, always be a cuboid and tri based geometries such
      // as a tetrahedron will be embedded within the cuboid of the array.
      // Hence, there are elements that will never be used and must be
      // marked accordingly
      std::fill(ijk2l_map_[var].begin(), ijk2l_map_[var].end(), -1);

      for (; !ijk_iter->is_at_end(); ++(*ijk_iter)) {
        l2ijk_map_[var][ijk_iter->get_current_local_id()] =
            ijk_iter->get_current_ijk_id();

        set_ijk2l_map_element(var, ijk_iter->get_current_ijk_id(),
                              ijk_iter->get_current_local_id());
      }

      assert(l2ijk_map_[var].size() == element->nb_dof(var));
    }
  }

  /// Translates a global dof id to a local dof id relative to the cell.
  /// @return whether the translation was successful
  /// @param global_id the global id that shall be translated
  /// @param local_id After a successful translation, contains the local
  /// id that corresponds to the supplied global id. Otherwise, equals -1.

  bool map_g2l(int global_id, int &local_id) const {
    UMAP_NAMESPACE_GMG::unordered_map< int, int >::const_iterator it =
        g2lmap_.find(global_id);

    if (it != g2lmap_.end()) {
      local_id = it->second;
      return true;
    }
    // The global id is invalid for this cell
    local_id = -1;
    return false;
  }

  /// Translates a vector of global dof ids to a vector with
  /// the corresponding local dof ids relative to the cell.
  /// @return true if the translation of all dofs was successful, or false
  /// if at least one translation has failed.
  /// @param global_ids A vector of global ids that shall be translated
  /// @param local_id After a successful translation, contains the local
  /// ids that corresponds to the supplied global ids. Entries for which
  /// the translation has failed will equal -1.

  bool map_g2l(const std::vector< int > &global_ids,
               std::vector< int > &local_ids) const {
    local_ids = std::vector< int >(global_ids.size(), -1);

    bool result = true;

    for (int i = 0; i < global_ids.size(); ++i) {
      // return false if at least one map_g2l call returned false
      if (result)
        result = map_g2l(global_ids[i], local_ids[i]);
      else
        map_g2l(global_ids[i], local_ids[i]);
    }

    return result;
  }

  /// Translates a local dof id to a global dof id.
  /// @return whether the translation was successful
  /// @param local_id the local id that shall be translated
  /// @param global_id After a successful translation, contains the global
  /// id that corresponds to the supplied local id. Otherwise, equals -1.

  bool map_l2g(int local_id, int var, int &global_id) const {
    global_id = element_->get_space().dof().cell2global(var, element_->cell_index(),		//TODO: cell2global correct?
                                               local_id);
    return true;
  }

  /// Translates a local dof id to a ijk dof id
  /// @return whether the translation was successful
  /// @param local_id The local dof id that shall be translated
  /// @param var The variable to which the local dof id belongs
  /// @param id After a successful translation, contains the ijk dof id
  /// that corresponds to the supplied local id.
  /// \c id[0] will contain the i-index, \c id[1] the j-index and \c id[2]
  /// the k index (if the cell is three dimensional). If the translation has
  /// failed, id will be empty.

  bool map_l2ijk(int local_id, int var, ijk_dof_id &id) const {
    assert(var >= 0 && var < l2ijk_map_.size());

    if (local_id < 0 || local_id >= l2ijk_map_[var].size()) {
      id.clear();
      return false;
    }

    id = l2ijk_map_[var][local_id];
    return true;
  }

  /// Maps a vector of local dof ids to their corresponding ijk ids.
  /// @return Whether the translation was successful. Returns false,
  /// if the translation of at least one dof id has failed.
  /// @param local_dof_ids A vector containing local dof ids that shall
  /// be translated.
  /// @param var The variable to which the local dof ids belong
  /// @param ids The i-th entry will contain the ijk id corresponding
  /// to the local dof id at \c local_dof_ids[i], or an empty ijk
  /// id if the translation has failed.

  bool map_l2ijk(const std::vector< int > &local_dof_ids, int var,
                 std::vector< ijk_dof_id > &ids) const {
    ids = std::vector< ijk_dof_id >(local_dof_ids.size());

    bool result = true;

    for (std::size_t i = 0; i < local_dof_ids.size(); ++i) {
      if (result)
        result = map_l2ijk(local_dof_ids[i], var, ids[i]);
      else
        map_l2ijk(local_dof_ids[i], var, ids[i]);
    }

    return result;
  }

  /// Translates an ijk dof id to its corresponding local dof id.
  /// @return whether the translation was successful
  /// @param id The ijk dof id that shall be translated
  /// @param var The variable to which the ijk dof id \c id belongs
  /// @param local_id Will contain the local dof id corresponding
  /// to the supplied ijk id, or -1 if the translation has failed.

  bool map_ijk2l(const ijk_dof_id &id, int var, int &local_id) const {
    assert(var >= 0);
    assert(id.size() == 2 || id.size() == 3);

    local_id = -1;

    if (id[0] >= ijk2l_map_[var].get_extent_of_dimension(0) || id[0] < 0)
      return false;
    if (id[1] >= ijk2l_map_[var].get_extent_of_dimension(1) || id[1] < 0)
      return false;
    if (id.size() == 3)
      if (id[2] >= ijk2l_map_[var].get_extent_of_dimension(2) || id[2] < 0)
        return false;

    local_id = get_ijk2l_map_element(var, id);

    if (local_id == -1)
      return false;
    return true;
  }

  /// Reinitializes the converter for a new cell. Only the data structures
  /// for translation between global and local elements depend on the actual
  /// cell object for which the \c DofIdConverter object was constructed.
  /// Thus, it is usually cheaper to use this function to rebind the converter
  /// to a new cell than to create a new DofIdConverter object if dof
  /// translation for a different cell is needed.
  /// @param element The new cell for which this \c DofIdConverter
  /// is supposed to translate dofs.
  /// \c element must be from the same mesh as the previously used cell
  /// (more specifically: must at least have the same cell type and the same
  /// fe degree for each variable)

  void rebind_to_element(const Element< DataType, DIM > *element) {
    assert(element != NULL);

    if (element->nb_fe() != this->element_->nb_fe()) //nb_fe correct?
      throw std::invalid_argument("Tried to rebind DofIdConverter to element"
                                  "with different number of variables");

    if (element->get_cell().cell_type().tag() != this->tag_)
      throw std::invalid_argument("Tried to rebind DofIdConverter to element"
                                  "of different type");

    for (int var = 0; var < this->element_->nb_fe(); ++var)

    {
      DofContainerLagrange <DataType, DIM > const * dof_element = dynamic_cast < DofContainerLagrange <DataType, DIM > const * >(element->get_fe(var)->dof_container());
          
      DofContainerLagrange <DataType, DIM > const * dof_element_ = dynamic_cast < DofContainerLagrange <DataType, DIM > const * >(element_->get_fe(var)->dof_container());
      assert(dof_element != nullptr);
      assert(dof_element_ != nullptr);
      
      if (dof_element->get_degree() !=
          dof_element_->get_degree())

        throw std::invalid_argument("Tried to rebind DofIdConverer to"
                                    "element with different fe degree");
    }

    element_ = element;

    init_g2l_map();
  }

  /// @return The cell type of the cell for which this converter was
  // constructed

  CellType::Tag get_cell_type() const { return tag_; }

  /// @return A pointer to the element for which this dof converter
  /// can convert ids

  const Element< DataType, DIM > *get_element() const { return element_; }

private:
  /// Initializes the global to local dof translation data structure

  void init_g2l_map() {
    int cell_index = element_->cell_index();

    g2lmap_.clear();

    for (int var = 0; var < element_->nb_fe(); ++var) {
      for (int local_dof = 0; local_dof < element_->nb_dof(var);
           ++local_dof) {

        int global_dof =
            element_->get_space().dof().cell2global(var, cell_index, local_dof); //TODO cell2global correct?

        g2lmap_.insert(std::make_pair(global_dof, local_dof));
      }
    }
  }

  /// Assigns a value to an entry in the ijk to local dof translation
  /// data structure. This function does not perform
  /// any bounds checks on the ijk id.
  /// @param var The variable to which the supplied local id belongs
  /// @param ijk_id The ijk dof id, of which the corresponding entry in
  /// the ijk to local dof translation data structure shall be modified
  /// @param local The local dof id, that shall be written to the
  /// appropriate position in the map. After a call to this function,
  /// \c ijk2l_map_[var][i][j][k]=local (in 3d) or \c
  /// ijk2l_map_[var][i][j]=local (in 2d)

  inline void set_ijk2l_map_element(int var, const ijk_dof_id &ijk_id,
                                    int local) {
    assert(ijk_id.size() == 2 || ijk_id.size() == 3);
    assert(ijk_id.size() == ijk2l_map_[var].get_dimension());
    assert(var >= 0 && var < ijk2l_map_.size());

    ijk2l_map_[var][ijk_id] = local;
  }

  /// Looks up the local id in the ijk to local dof translation data structure
  /// that corresponds to a given ijk id. This function does not perform
  /// any bounds checks on the ijk id.
  /// @param var The variable to which the supplied ijk dof id belongs
  /// @param ijk_id The ijk id that shall be looked up
  /// @return The local id that corresponds to the given ijk id

  inline int get_ijk2l_map_element(int var, const ijk_dof_id &ijk_id) const {
    assert(ijk_id.size() == 2 || ijk_id.size() == 3);
    assert(ijk_id.size() == ijk2l_map_[var].get_dimension());
    assert(var >= 0 && var < ijk2l_map_.size());

    return ijk2l_map_[var][ijk_id];
  }

  const Element< DataType, DIM > *element_;
  CellType::Tag tag_;
};

/// This class creates and caches \c DofIdConverter objects. When the user
/// requests access to a \c DofIdConverter for a specific cell type,
/// this class will check if it has created a converter for this cell type
/// before and, if so, will rebind this converter to the cell for which
/// a converter has been requested. Otherwise a new converter will be created.
/// This is a useful performance optimization when iterating through many
/// cells of the same type because it will avoid recreation and full constructor
/// initialization of dof converter objects.

template < class LAD, int DIM > 
class OnDemandDofIdConverter {
public:
  TYPE_FROM_CLASS(LAD, DataType);

  ~OnDemandDofIdConverter() {
    for (int i = 0; i < converters_.size(); ++i)
      if (converters_[i] != NULL)
        delete converters_[i];
  }

  /// Obtain a dof id converter for a given cell. When a dof id converter
  /// has been returned before for a cell of the same type, the previously
  /// created converter will be rebound to the new cell, otherwise a
  /// new converter will be created.
  /// @param element The cell for which a converter is required. Unless you
  /// really know what you are doing, all cells that are used as a parameter
  /// here must originate from the same mesh.
  /// @return A reference to a \c DofIdConverter for the given cell.

  DofIdConverter< LAD, DIM > &get_converter(const Element< DataType, DIM > *element) {
    CellType::Tag t = element->get_cell().cell_type().tag();

    for (int i = 0; i < converters_.size(); ++i) {
      if (converters_[i] != NULL) {
        if (converters_[i]->get_cell_type() == t) {
          converters_[i]->rebind_to_element(element);
          return *(converters_[i]);
        }
      }
    }
    converters_.push_back(new DofIdConverter< LAD, DIM >(element));
    return *(converters_.back());
  }

private:
  std::vector< DofIdConverter< LAD, DIM > * > converters_;
};

/// This class (and its derived classes) abstract the lookup if a dof is unknown

template < class LAD, int DIM > 
class IsDofKnownLookup {
public:
  /// construct object
  /// @param dof_ident A shared pointer to a DofIdentification object that
  /// contains information about which dofs correspond to a dof on a coarser
  /// mesh.

  explicit IsDofKnownLookup(
      const boost::shared_ptr< const DofIdentification< LAD, DIM > > &dof_ident)
      : dof_ident_(dof_ident) {}

  /// @return whether a dof is known
  /// @param dof_id The global id of the dof that shall be looked up
  virtual bool operator()(const int dof_id) const = 0;

protected:
  boost::shared_ptr< const DofIdentification< LAD, DIM > > dof_ident_;
};

/// Implements a direct lookup via the DofIdentification. This means, that
/// DoFs are considered known iff the DofIdentification has identified
/// them with dofs on a coarser grid.

template < class LAD, int DIM >
class DirectIsDofKnownLookup : public IsDofKnownLookup< LAD, DIM > {
public:
  /// construct object
  /// @param dof_ident A shared pointer to a DofIdentification object that
  /// contains information about which dofs correspond to a dof on a coarser
  /// mesh.

  explicit DirectIsDofKnownLookup(
      const boost::shared_ptr< const DofIdentification< LAD, DIM > > &dof_ident)
      : IsDofKnownLookup< LAD, DIM >(dof_ident) {}

  /// @return whether a dof is known
  /// @param dof_id The global id of the dof that shall be looked up

  virtual bool operator()(const int dof_id) const {
    return this->dof_ident_->dof_exists_on_coarse_level(dof_id);
  }
};

/// Implements a lookup via an overlay. This means, that these objects contain
/// a bitmap in which dofs can marked as known. A dof will then be considered
/// known, if it is either marked as known in the overlay bitmap or considered
/// known by the DofIdentification.

template < class LAD, int DIM >
class OverlayBasedIsDofKnownLookup : public IsDofKnownLookup< LAD, DIM > {
public:
  TYPE_FROM_CLASS(LAD, DataType);
  /// construct object
  /// @param dof_ident A shared pointer to a DofIdentification object that
  /// contains information about which dofs correspond to a dof on a coarser
  /// mesh.
  /// @param dof_converter A \c DofIdConverter object for the cell on which
  /// this \c OverlayBasedIsDofKnownLookup object will be used

  OverlayBasedIsDofKnownLookup(
      const boost::shared_ptr< const DofIdentification< LAD, DIM > > &dof_ident,
      const DofIdConverter< LAD, DIM > *dof_converter)
      : IsDofKnownLookup< LAD, DIM >(dof_ident), dof_converter_(dof_converter),
        use_overlay_(false) {}

  /// @return whether a dof is known
  /// @param dof_id The global id of the dof that shall be looked up

  virtual bool operator()(const int dof_id) const {
    int local_id = -1;
    if (use_overlay_) {
      if (dof_converter_->map_g2l(dof_id, local_id)) {
        assert(local_id >= 0 && local_id < overlay_.size());

        if (overlay_[local_id])
          return true;
      }
    }
    return this->dof_ident_->dof_exists_on_coarse_level(dof_id);
  }

  /// marks a dof with a given local id (with respect to the cell) as known
  /// @param local_id The local dof id that shall be marked as known

  void mark_dof_as_known(int local_id) {
    assert(local_id >= 0 && local_id < overlay_.size());
    overlay_[local_id] = true;
  }

  /// unmarks a dof with a given local id (with respect to the cell) as known
  /// @param local_id The local dof id that shall be unmarked as known

  void unmark_dof_as_known(int local_id) {
    assert(local_id >= 0 && local_id < overlay_.size());
    overlay_[local_id] = false;
  }

  /// Sets the size of the overlay to the number of dofs of a specified
  /// variable. When the overlay shall be used, a call to this method is
  /// neccessary because it creates the needed space for a bitmap for dofs of
  /// the specified variable.
  /// @param element The element on which the overlay shall be used
  /// @param var The variable id for which the overlay shall be created

  void resize_overlay_to_var(const Element< DataType, DIM > &element, int var) {
    overlay_ = std::vector< bool >(element.nb_dof(var), false);			//TODO var=fe_ind?
  }

  /// Enables the overlay. For the overlay to function, a preceding call
  /// to \c resize_overlay_to_var is required.

  void enable_overlay() { use_overlay_ = true; }

  /// Disables the overlay

  void disable_overlay() { use_overlay_ = false; }

  /// Set if the overlay shall be used
  /// @param enabled whether the overlay shall be used.

  void use_overlay(bool enabled) { use_overlay_ = enabled; }

private:
  bool use_overlay_;
  const DofIdConverter< LAD, DIM > *dof_converter_;
  std::vector< bool > overlay_;
};

/// Finds the nearest known dofs relative to a given dof in a cell

template < class LAD, int DIM > 
class NearestKnownDofsFinder {
public:
  TYPE_FROM_CLASS(LAD, VectorType);
  TYPE_FROM_CLASS(LAD, DataType);

  typedef Vec<DIM, typename LAD::DataType> Coord;

  /// @param dof_converter A pointer to a \c DofIdConverter object
  /// for the cell that has been supplied in \c element
  /// @param lookup A pointer to a \c IsDofKnownLookup that will be
  /// used to check if a Dof is known

  explicit NearestKnownDofsFinder(const DofIdConverter< LAD, DIM > *dof_converter,
                                  const IsDofKnownLookup< LAD, DIM > *lookup)
      : dof_converter_(dof_converter), is_dof_known_(lookup) {}

  /// finds the nearest known dofs around a given a local dof, and returns
  /// a vector containing the global dof ids of these nearest dofs.
  /// @param local_dof_id The local dof id (relative to the cell and variable)
  /// around which nearest known dofs shall be searched
  /// @param var The variable to which \c local_dof_id belongs
  /// @param nearest_dofs A vector that will contain the global dof ids
  /// of the nearest dof ids that have been found once this function call has
  /// returned.

  void find_nearest_known_dofs(int local_dof_id, int var,
                               std::vector< int > &nearest_dofs) const;

  /// Finds the nearest known dof on a subentity
  /// @return The global dof id of the nearest known dof or -1, if no
  /// known dof has been found.
  /// @param tdim The dimension of the subentity
  /// @param sindex The sindex of the subentity (unused if tdim == tdim of cell)
  /// @param global_dof The global dof of which the nearest dofs shall be found.
  /// It must also reside on the supplied subentity.
/*
  int find_nearest_known_dof_on_subent(TDim tdim, int sindex, int global_dof,
                                       int var) const;
*/
  /// contains the ijk variations that are used to find nearest dofs on Q cells

  struct QDeltas {
    static const int deltas2d_on_axis[4][3];

    static const int deltas2d_diagonals[4][3];

    static const int deltas3d_on_axis[6][3];

    static const int deltas3d_two_axis_diagonals[12][3];

    static const int deltas3d_three_axis_diagonals[8][3];
  };

  /// contaons the ijk variations that are used to find nearest dofs on P cells

  struct PDeltas {
    static const int deltas2d[6][3];

    static const int deltas3d[12][3];
  };

private:
  /// Processes a list of ijk variations by applying them to a ijk dof id
  /// representing the origin of the search, checking if the dofs are known,
  /// and if so, inserting them into a result vector
  /// @param delta_list The ijk variation list
  /// @param origin A ijk dof id representing the origin of the search
  /// @param var The variable to which the ijk dof id \c origin belongs
  /// @param found_dofs A vector that will contain the nearest known dofs
  /// that have been found once a call to this function returns.

  template < unsigned N >
  inline void process_delta_list(const int (&delta_list)[N][3],
                                 ijk_dof_id &origin, int var,
                                 std::vector< int > &found_dofs) const {
    std::size_t num_deltas = N;
    for (std::size_t i = 0; i < num_deltas; ++i)
      investigate_at_delta(origin, delta_list[i][0], delta_list[i][1],
                           delta_list[i][2], var, found_dofs);
  }

  /// Performs the search for nearest known dofs on Q (quad based) cells
  /// @param origin A ijk dof id representing the origin of the search
  /// @param var The variable to which the ijk dof id \c origin belongs
  /// @param found_dofs A vector that will contain the nearest known dofs
  /// that have been found once a call to this function returns.

  void q_proximity_search(ijk_dof_id &search_origin, int var,
                          std::vector< int > &nearest_dofs) const;

  /// Performs the search for nearest known dofs on P (tri based) cells
  /// @param origin A ijk dof id representing the origin of the search
  /// @param var The variable to which the ijk dof id \c origin belongs
  /// @param found_dofs A vector that will contain the nearest known dofs
  /// that have been found once a call to this function returns.

  void p_proximity_search(ijk_dof_id &search_origin, int var,
                          std::vector< int > &nearest_dofs) const {
    assert(search_origin.size() == 2 || search_origin.size() == 3);

    if (search_origin.size() == 2)
      process_delta_list(PDeltas::deltas2d, search_origin, var, nearest_dofs);
    else
      process_delta_list(PDeltas::deltas3d, search_origin, var, nearest_dofs);
  }

  /// Applies an ijk variation ("delta") and checks if the resultung dof
  /// is known. If so, it is inserted into a vector.
  /// @param origin An ijk dof id on which the variation shall be applied
  /// @param delta_i The variation in i direction
  /// @param delta_j The variation in j direction
  /// @param delta_k The variation in k direction
  /// @param var The variable to which \c origin belongs
  /// @param found_dofs A vector into which this function will insert
  /// the global id of the dof that results from the variation if it is known.

  inline void investigate_at_delta(ijk_dof_id &origin, int delta_i, int delta_j,
                                   int delta_k, int var,
                                   std::vector< int > &found_dofs) const {
    assert(origin.size() == 2 || origin.size() == 3);

    if (origin.size() == 2)
      assert(delta_k == 0);

    origin[0] += delta_i;
    origin[1] += delta_j;
    origin[2] += delta_k;

    investigate_dof(origin, var, found_dofs);

    origin[0] -= delta_i;
    origin[1] -= delta_j;
    origin[2] -= delta_k;
  }

  /// Checks if a given dof id is known and if so, stores the corresponding
  /// global id in a std::vector
  /// @param dof_id The ijk dof id that shall be checked
  /// @param var The variable to which the ijk id belongs
  /// @param found_dofs A vector into which this function will insert
  /// the global id of the dof if it is known.

  inline void investigate_dof(const ijk_dof_id &dof_id, int var,
                              std::vector< int > &found_dofs) const {
    int global_id = -1;
    if (is_dof_known(dof_id, var, global_id))
      found_dofs.push_back(global_id);
  }

  /// Checks if an ijk dof is known, and stores the global id that
  /// is associated with the ijk dof
  /// @return whether the given dof is known, ie also exists on the next
  /// coarser grid
  /// @param dof_id The ijk dof that shall be checked
  /// @param var The variable to which the supplied ijk dof belongs
  /// @param found_global_id After a call to this function, will either
  /// equal the global id that corresponds to the supplied ijk dof id,
  /// or will equal -1 if the translation from ijk to global id has failed.

  inline bool is_dof_known(const ijk_dof_id &dof_id, int var,
                           int &found_global_id) const {
    found_global_id = -1;
    int local_id = -1;
    if (!dof_converter_->map_ijk2l(dof_id, var, local_id))
      return false;

    if (!dof_converter_->map_l2g(local_id, var, found_global_id))
      return false;

    return (*is_dof_known_)(found_global_id);
  }

  const DofIdConverter< LAD, DIM > *dof_converter_;
  const IsDofKnownLookup< LAD, DIM > *is_dof_known_;
};

template < class LAD, int DIM >
const int NearestKnownDofsFinder< LAD, DIM >::QDeltas::deltas2d_on_axis[4][3] = {
    // on axis neighbors
    {+1, 0, 0},
    {-1, 0, 0},
    {0, +1, 0},
    {0, -1, 0}};

template < class LAD, int DIM >
const int NearestKnownDofsFinder< LAD, DIM >::QDeltas::deltas2d_diagonals[4][3] = {
    // diagonals
    {+1, +1, 0},
    {-1, -1, 0},
    {-1, +1, 0},
    {+1, -1, 0}};

template < class LAD, int DIM >
const int NearestKnownDofsFinder< LAD, DIM >::QDeltas::deltas3d_on_axis[6][3] = {
    // on axis neighbors
    {+1, 0, 0}, {-1, 0, 0}, {0, +1, 0}, {0, -1, 0}, {0, 0, +1}, {0, 0, -1}};

template < class LAD, int DIM >
const int NearestKnownDofsFinder<
    LAD, DIM >::QDeltas::deltas3d_two_axis_diagonals[12][3] = {
    // two axis diagonals
    {+1, +1, 0}, {-1, -1, 0}, {-1, +1, 0}, {+1, -1, 0},
    {+1, 0, +1}, {-1, 0, -1}, {-1, 0, +1}, {+1, 0, -1},
    {0, +1, +1}, {0, -1, -1}, {0, -1, +1}, {0, +1, -1}};

template < class LAD, int DIM >
const int NearestKnownDofsFinder<
    LAD, DIM >::QDeltas::deltas3d_three_axis_diagonals[8][3] = {
    // three axis diagonals
    {+1, +1, +1}, {-1, -1, +1}, {-1, +1, +1}, {+1, -1, +1},
    {+1, +1, -1}, {-1, -1, -1}, {-1, +1, -1}, {+1, -1, -1}};

template < class LAD, int DIM >
const int NearestKnownDofsFinder< LAD, DIM >::PDeltas::deltas2d[6][3] = {
    // on axis neighbors
    {+1, 0, 0},
    {-1, 0, 0},
    {0, +1, 0},
    {0, -1, 0},
    // antisymmetric diagonals
    {+1, -1, 0},
    {-1, +1, 0}};

template < class LAD, int DIM >
const int NearestKnownDofsFinder< LAD, DIM >::PDeltas::deltas3d[12][3] = {
    // on axis neighbors
    {+1, 0, 0},
    {-1, 0, 0},
    {0, +1, 0},
    {0, -1, 0},
    // antisymmetric diagonals
    {+1, -1, 0},
    {-1, +1, 0},

    // k =/= 0
    // on axis neighbors
    {0, 0, +1},
    {0, 0, -1},
    // antisymmetric diagonals
    {+1, 0, -1},
    {-1, 0, +1},
    {0, +1, -1},
    {0, -1, +1}};

template < class LAD, int DIM >
void NearestKnownDofsFinder< LAD, DIM >::find_nearest_known_dofs(
    int local_dof_id, int var, std::vector< int > &nearest_dofs) const {
  nearest_dofs.clear();

  const Element< DataType, DIM > *element = dof_converter_->get_element();

  ijk_dof_id search_origin;
  if (dof_converter_->map_l2ijk(local_dof_id, var, search_origin)) {
    if (element->get_cell().cell_type().tag() == CellType::QUADRILATERAL ||
        element->get_cell().cell_type().tag() == CellType::HEXAHEDRON)
      q_proximity_search(search_origin, var, nearest_dofs);
    else
      p_proximity_search(search_origin, var, nearest_dofs);
  }
}

/// Finds the nearest known dof on a subentity
/// @return The global dof id of the nearest known dof or -1, if no
/// known dof has been found.
/// @param tdim The dimension of the subentity
/// @param sindex The sindex of the subentity (unused if tdim == tdim of cell)
/// @param global_dof The global dof of which the nearest dofs shall be found.
/// It must also reside on the supplied subentity.
 /*
template < class LAD, int DIM >
int NearestKnownDofsFinder< LAD, DIM >::find_nearest_known_dof_on_subent(
    TDim tdim, int sindex, int global_dof, int var) const {
  std::vector< Coord > coordinates;
  std::vector< int > dofs;

  const Element< DataType, DIM > *element = dof_converter_->get_element();

  DofContainerLagrange <DataType, DIM > const * dof_element = 
  dynamic_cast < DofContainerLagrange <DataType, DIM > const * >(element->get_fe(var)->dof_container());




  if (tdim < element->get_cell().tdim()) {
    coordinates = dof_element.get_dof_coords_on_subentity(tdim, sindex);
    
    dofs = dof_element.get_dof_on_subentity(tdim, sindex);
    
  } else {
    dofs = dof_element->get_dof_coords
    element->aaa_get_dof_indices(var, dofs);
    element->space().dof().aaa_get_coord_on_cell(var, element->get_cell_index(),
                                             coordinates);
  }

  double min_distance = std::numeric_limits< double >::max();
  int min_dof = -1;

  assert(coordinates.size() == dofs.size());

  Coord search_dof_coord;
  bool found = false;
  for (std::size_t i = 0; i < dofs.size() && (!found); ++i) {
    if (dofs[i] == global_dof)
    {
      search_dof_coord = coordinates[i];
      found = true;
    }
  }

  if (!found)
    // We have not found the dof we are looking for in the coordinate list
    return -1;

  for (std::size_t i = 0; i < coordinates.size(); ++i) {
    assert(coordinates[i].size() == search_dof_coord.size());

    if (dofs[i] != global_dof && (*is_dof_known_)(dofs[i])) {
      double distance = 0.0;

      // calculate distance
      for (std::size_t axis = 0; axis < search_dof_coord.size(); ++axis)
        distance += std::abs(search_dof_coord[axis] - coordinates[i][axis]);

      if (distance < min_distance) {
        min_distance = distance;
        min_dof = dofs[i];
      }
    }
  }

  return min_dof;
}*/

/// Performs the search for nearest known dofs on Q (quad based) cells
/// @param origin A ijk dof id representing the origin of the search
/// @param var The variable to which the ijk dof id \c origin belongs
/// @param found_dofs A vector that will contain the nearest known dofs
/// that have been found once a call to this function returns.

template < class LAD, int DIM >
void NearestKnownDofsFinder< LAD, DIM >::q_proximity_search(
    ijk_dof_id &search_origin, int var,
    std::vector< int > &nearest_dofs) const {
  assert(search_origin.size() == 2 || search_origin.size() == 3);

  if (search_origin.size() == 2) {
    // 2d case

    // on axis neighbors
    process_delta_list(QDeltas::deltas2d_on_axis, search_origin, var,
                       nearest_dofs);

    if (nearest_dofs.empty())
      // look at diagonals
      process_delta_list(QDeltas::deltas2d_diagonals, search_origin, var,
                         nearest_dofs);
  } else {
    // 3d case

    // First, try the on axis neighbors
    process_delta_list(QDeltas::deltas3d_on_axis, search_origin, var,
                       nearest_dofs);

    if (nearest_dofs.empty()) {
      // Then look at two axis diagonals
      process_delta_list(QDeltas::deltas3d_two_axis_diagonals, search_origin,
                         var, nearest_dofs);

      if (nearest_dofs.empty())
        // Finally, look at three axis diagonals
        process_delta_list(QDeltas::deltas3d_three_axis_diagonals,
                           search_origin, var, nearest_dofs);
    }
  }
}

/// Iterates through all cells of a mesh and executes a function for
/// each cell. The iteration is guaranteed to be in ascending order with
/// respect to the cell indices.

template < class LAD, int DIM > 
class ForEachCell {
public:
  TYPE_FROM_CLASS(LAD, MatrixType);
  TYPE_FROM_CLASS(LAD, DataType);

  typedef boost::shared_ptr< VectorSpace< DataType, DIM > > SpacePtr;
  typedef boost::shared_ptr< MatrixType > MatrixPtr;

  /// Construct object
  /// @param mesh The mesh to iterate
  /// @param space The vectorspace on which the finite element functions
  /// live

  ForEachCell(MeshPtr mesh, SpacePtr space) : mesh_(mesh), space_(space) {
    assert(mesh_ != NULL);
    assert(space_ != NULL);
  }

  /// Execute the iteration.
  /// @param f The function that is executed for each cell. It is required
  /// to accept one Parameter of type \c Element<Scalar>& or \c const
  /// Element<Scalar>&.

  template < class Function > void operator()(Function &f) {
    assert(mesh_ != NULL);
    assert(space_ != NULL);

    unsigned num_cells = mesh_->num_entities(mesh_->tdim());

    for (unsigned cell_index = 0; cell_index < num_cells; ++cell_index) {
      Element< DataType, DIM > current_cell(*space_, cell_index);
      f(current_cell);
    }
  }

  template < class Function >
  void operator()(Function &f, MatrixPtr matrix, std::set< int > &set) {
    assert(mesh_ != NULL);
    assert(space_ != NULL);
    assert(matrix != NULL);

    unsigned num_cells = mesh_->num_entities(mesh_->tdim());

    for (unsigned cell_index = 0; cell_index < num_cells; ++cell_index) {
      Element< DataType, DIM > current_cell(*space_, cell_index);
      f(current_cell, matrix.get(), set);
    }
  }

private:
  MeshPtr mesh_;
  SpacePtr space_;
};

/**  THIS CODE IS EXPERIMENTAL AND SHOULD NOT BE USED IN ANY REAL APPLICATION
class BestSurroundingDofSubsetCalculator
{
  typedef UMAP_NAMESPACE_GMG::unordered_map<int, Coord> DofCoordMapType;
public:

  BestSurroundingDofSubsetCalculator()
  {
  }

  BestSurroundingDofSubsetCalculator(const Element<DataType>& cell)
  {
    load_dofs(cell);
  }

  void load_dofs(const Element<DataType>& cell)
  {
    int cell_tdim = cell.get_cell().tdim();
    int num_cells = cell.get_cell().num_incident_entities(cell_tdim);
    int num_vars = cell.get_num_variables();
    int total_num_dofs = 0;
    for (int var = 0; var < num_vars; ++var)
      total_num_dofs += cell.get_num_dofs(var);

    dof_coord_map_ = UMAP_NAMESPACE_GMG::unordered_map<int, Coord>(1.2 *
total_num_dofs * (num_cells + 1));

    for (int var = 0; var < num_vars; ++var)
    {
      for (IncidentEntityIterator neighbor_cell =
cell.get_cell().begin_incident(cell_tdim); neighbor_cell !=
cell.get_cell().end_incident(cell_tdim); ++neighbor_cell)
      {
        int cell_idx = neighbor_cell->index();
        load_dofs_of_cell(cell.space().dof(), var, cell_idx);
      }

      load_dofs_of_cell(cell.space().dof(), var, cell.get_cell_index());
    }
  }

  void operator()(const std::vector<int>& dofs,
          int center_dof_global_id,
          std::vector<int>& out) const
  {
    std::vector<Coord> coordinates;
    coordinates.reserve(dofs.size());
    for (int i = 0; i < dofs.size(); ++i)
    {
      DofCoordMapType::const_iterator it = dof_coord_map_.find(dofs[i]);
      assert(it != dof_coord_map_.end());

      coordinates.push_back(it->second);
    }

    DofCoordMapType::const_iterator coords_of_center_dof =
dof_coord_map_.find(center_dof_global_id); assert(coords_of_center_dof !=
dof_coord_map_.end());

    (*this)(dofs, coordinates, coords_of_center_dof->second, out);
  }

  void operator()(const std::vector<int>& dofs,
          const std::vector<Coord>& coordinates,
          const Coord& desired_barycenter,
          std::vector<int>& out) const
  {
    assert(dofs.size() == coordinates.size());

    out.clear();

    if (dofs.size() == 0)
      return;

    out.reserve(coordinates.size());

    // Specifies which dofs are used for the barycenter calculation
    std::vector<bool> used_points_map(coordinates.size(), true);

    // optimize as long as possible
    while (optimize_used_map(coordinates, desired_barycenter, used_points_map))
      ;

    for (std::size_t i = 0; i < used_points_map.size(); ++i)
      if (used_points_map[i])
        out.push_back(dofs[i]);
#ifdef DEBUG
    std::cout << "new #dofs: " << out.size() << " old: " << dofs.size() <<
std::endl; #endif
  }

private:

  void load_dofs_of_cell(const DofPartition<DataType>& dof_partition, int var,
int cell_idx)
  {
    std::vector<int> dofs_on_cell;
    std::vector<Coord> coords;

    dof_partition.get_coord_on_cell(var, cell_idx, coords);
    dof_partition.get_dofs_on_cell(var, cell_idx, dofs_on_cell);

    assert(dofs_on_cell.size() == coords.size());

    for (std::size_t i = 0; i < dofs_on_cell.size(); ++i)
      dof_coord_map_.insert(std::make_pair(dofs_on_cell[i], coords[i]));
  }

  inline bool optimize_used_map(const std::vector<Coord>& coordinates,
                                const Coord& desired_barycenter,
                                std::vector<bool>& used_map) const
  {
    if (coordinates.size() == 0)
      return false;

    Coord calculated_barycenter;
    calculate_barycenter(coordinates, used_map, calculated_barycenter);

    double min_distance = distance(calculated_barycenter, desired_barycenter);
    int min_index = -1;

    for (std::size_t i = 0; i < used_map.size(); ++i)
    {
      // if we use this dof for barycenter calculation try disabling
      // it and see if we find a new minimum
      if (used_map[i])
      {
        used_map[i] = false;

        calculate_barycenter(coordinates, used_map, calculated_barycenter);
        double current_distance = distance(calculated_barycenter,
desired_barycenter);

        if (current_distance < min_distance)
        {
          min_distance = current_distance;
          min_index = i;
        }

        used_map[i] = true;
      }
    }

    if (min_index != -1)
    {
      used_map[min_index] = false;
      return true;
    }
    return false;
  }

  void calculate_barycenter(const std::vector<Coord>& coords,
                            const std::vector<bool>& used_points,
                            Coord& barycenter) const
  {
    assert(coords.size() == used_points.size());

    if (coords.size() == 0)
    {
      barycenter = Coord();
      return;
    }

    std::size_t dimensions = coords.front().size();

    barycenter = Coord(dimensions, 0.0);
    std::size_t num_contributions = 0;

    for (std::size_t i = 0; i < coords.size(); ++i)
    {
      if (used_points[i])
      {
        for (std::size_t j = 0; j < dimensions; ++j)
          barycenter[j] += coords[i][j];
        ++num_contributions;
      }
    }

    if (num_contributions != 0)
      for (std::size_t j = 0; j < dimensions; ++j)
        barycenter[j] /= static_cast<double> (num_contributions);
  }

  inline double distance(const Coord& point0, const Coord& point1) const
  {
    assert(point0.size() == point1.size());
    double result = 0.0;
    for (std::size_t i = 0; i < point0.size(); ++i)
      result += std::abs(point0[i] - point1[i]);

    return result;
  }

  DofCoordMapType dof_coord_map_;
};
 END EXPERIMENTAL CODE */

/// Interpolates the unkown dofs on a cell.
/// This only works if the mesh on which the interpolation takes place is
/// not the coarsest mesh in the hierarchy. The algorithm used here can only
/// interpolate vectors that have been transferred from the next coarser grid.

template < class LAD, int DIM > 
class LinearCellInterpolator {
public:
  TYPE_FROM_CLASS(LAD, VectorType);
  TYPE_FROM_CLASS(LAD, DataType);

  /// Construct the interpolator.
  /// @param dof_ident A pointer to a dof identification object that identifies
  /// the dofs on the mesh on which the interpolation happens with the dofs on
  /// the next coarser mesh
  /// @param vector A pointer to the vector that shall be interpolated.
  /// This pointer will only be used to read the values of the known dofs.
  /// @param still_unknown_dofs A pointer to a bitmap where dofs that could
  /// not be interpolated will be flagged. If NULL, then the bitmap will be
  /// ignored.
  /// @param interpolate_from_interpolated If true, this cell interpolator
  /// will only attempt to interpolate dof that have been flagged as still
  /// unknown in \c still_unknown_dofs. The interpolation will then
  /// interpolate these dofs from interpolated dofs of a prior interpolation.

  LinearCellInterpolator(
      const boost::shared_ptr< const DofIdentification< LAD, DIM > > &dof_ident,
      VectorType *vector, std::vector< bool > *still_unknown_dofs,
      bool interpolate_from_interpolated = false)
      : dof_ident_(dof_ident), vector_(vector),
        still_unknown_map_(still_unknown_dofs),
        interpolate_from_interpolated_(interpolate_from_interpolated) {
    assert(vector != NULL);
  }

  /// Interpolates the unkown dofs on a cell.
  /// This only works if the mesh on which the interpolation takes place is
  /// not the coarsest mesh in the hierarchy. The algorithm used here can only
  /// interpolate vectors that have been transferred from the next coarser grid.

  void operator()(const Element< DataType, DIM > &element)
  {

      int cell_dimension = element.get_cell().tdim();

      const DofPartition< DataType, DIM > &dof_partition = element.get_space().dof();

      const DofIdConverter< LAD, DIM > &dof_converter =
          on_demand_dof_converter_.get_converter(&element);

      OverlayBasedIsDofKnownLookup< LAD, DIM > is_dof_known(dof_ident_, &dof_converter);

      // BestSurroundingDofSubsetCalculator subset_calculator(element);

      for (int var = 0; var < element.nb_fe(); ++var) {		//TODO: replace var by fe_ind?
        // allows us to modify the overlay
        is_dof_known.resize_overlay_to_var(element, var);
        // disable overlay by default
        is_dof_known.disable_overlay();

        if (interpolate_from_interpolated_) {
          // configure IsDofKnownLookup such that all dofs
          // are known except the ones that are marked in
          // the still unknown bitmap
          for (int local_dof_id = 0; local_dof_id < element.nb_dof(var);
               ++local_dof_id) {
            int global_id = -1;
            dof_converter.map_l2g(local_dof_id, var, global_id);

            if (dof_partition.is_dof_on_subdom(global_id)) {
              if (!is_local_dof_still_unknown(local_dof_id, var, element,
                                              dof_partition))
                is_dof_known.mark_dof_as_known(local_dof_id);
            } else
              is_dof_known.mark_dof_as_known(local_dof_id);
          }
        }

        std::vector< int > dofs;
        for (int dim = 0; dim < cell_dimension; ++dim) {
          int sindex = 0;
          for (IncidentEntityIterator ent = element.get_cell().begin_incident(dim);
               ent != element.get_cell().end_incident(dim); ++ent, ++sindex) {
            // obtain dofs
            dof_partition.get_dofs_on_subentity(var, element.cell_index(), dim,
                                                sindex, dofs);
            interpolate_entity(var, element, *ent, dofs, dof_converter,
                               is_dof_known);
          }
        }
        // for dim == cell_dimension
        element.get_dof_indices(var, dofs);
        interpolate_entity(var, element, element.get_cell(), dofs, dof_converter,
                           is_dof_known);
      }
    }

private:
  /// @return whether the supplied entity is adjacent to cells that have
  /// already been fully interpolated.
  /// @param ent The entity that shall be checked
  /// @param element The cell to which the entity belongs

  bool has_entity_interpolated_neighbor_cells(
      const Entity &ent, const Element< DataType, DIM > &element) const {
    int cell_tdim = element.get_cell().tdim();
    for (IncidentEntityIterator cell = ent.begin_incident(cell_tdim);
         cell != ent.end_incident(cell_tdim); ++cell)
      if (is_cell_already_interpolated(*cell, element.get_cell()))
        return true;

    return false;
  }

  /// Interpolate an entity
  /// @param var The variable id
  /// @param element The cell on which the entity resides
  /// @param ent The entity that shall be interpolated
  /// @param dofs The dofs that are on this entity
  /// @param dof_converter A DofIdConverter for the cell on which the entity
  /// resides
  /// @param is_dof_known An \c OverlayBasedIsDofKnownLookup that will be used
  /// to query if a dof is unknown, ie needs to be interpolated.

  void interpolate_entity(int var, const Element< DataType, DIM > &element,
                          const Entity &ent, const std::vector< int > &dofs,
                          const DofIdConverter< LAD, DIM > &dof_converter,
                          OverlayBasedIsDofKnownLookup< LAD, DIM > &is_dof_known);

  /// @return whether a global dof is still unknown according to the
  /// \c still_unknown_map. This corresponds to a lookup in the \c
  /// still_unknown_map bitmap. This only works if a non-NULL \c
  /// still_unknown_dofs vector has been supplied to the \c
  /// LinearCellInterpolator constructor.
  /// @param global_dof The global id of the dof that shall be looked up
  /// @param dof_partition A dof partition object to convert global dof ids
  /// to local ids (with respect to the subdomain)

  inline bool is_global_dof_still_unknown(
      int global_dof, const DofPartition< DataType, DIM > &dof_partition) const {
    assert(still_unknown_map_ != 0);

    int id_on_sd = 0;
    dof_partition.global2local(global_dof, &id_on_sd);
    if (still_unknown_map_ != NULL)
      return (*still_unknown_map_)[id_on_sd];
    else
      return false;
  }

  /// @return whether a local (with respect to the cell) dof is still unknown
  /// according to the \c still_unknown_map. This corresponds to a lookup in the
  /// \c still_unknown_map bitmap. This only works if a non-NULL \c
  /// still_unknown_dofs vector has been supplied to the \c
  /// LinearCellInterpolator constructor.
  /// @param local_dof The local id of the dof that shall be looked up
  /// @param var The variable to which the local dof id belongs
  /// @param cell The cell on which the local dof resides
  /// @param dof_partition A dof partition object to convert global dof ids
  /// to local ids (with respect to the subdomain)

  inline bool is_local_dof_still_unknown(
      int local_dof, int var, const Element< DataType, DIM > &cell,
      const DofPartition< DataType, DIM > &dof_partition) const {
    assert(still_unknown_map_ != 0);

    int global_id = dof_partition.cell2global(var, cell.cell_index(), local_dof);		//TODO: cell2global correct?
    return is_global_dof_still_unknown(global_id, dof_partition);
  }

  /// @return whether a given cell has already been fully interpolated
  /// @param cell The cell that shall be checked
  /// @param currently_interpolated_cell The cell that is currently intpolated

  inline bool is_cell_already_interpolated(
      const Entity &cell, const Entity &currently_interpolated_cell) const {
    return cell.index() < currently_interpolated_cell.index();
  }

  /// Sets the status of a dof in the \c still_unknown_map_.
  /// @param dof_partition A dof partition object used to map global dof ids
  /// to local (with respect to the subdomain) ids
  /// @param global_dof_id The global id of the dof of which the status shall
  /// be set
  /// @param still_unknown The value to which the dof's entry in the map will
  /// be set. true means that the dof could not be interpolated and is still
  /// unknown

  inline void
  mark_dof_in_still_unknown_map(const DofPartition< DataType, DIM > &dof_partition,
                                int global_dof_id, bool still_unknown) {
    assert(dof_partition.is_dof_on_subdom(global_dof_id));

    if (still_unknown_map_ != NULL) {
      int id_on_sd = 0;
      dof_partition.global2local(global_dof_id, &id_on_sd);

      (*still_unknown_map_)[id_on_sd] = still_unknown;
    }
  }

  bool interpolate_from_interpolated_;
  boost::shared_ptr< const DofIdentification< LAD, DIM > > dof_ident_;
  VectorType *vector_;
  OnDemandDofIdConverter< LAD, DIM > on_demand_dof_converter_;
  std::vector< bool > *still_unknown_map_;
};

/// Interpolate an entity
/// @param var The variable id
/// @param element The cell on which the entity resides
/// @param ent The entity that shall be interpolated
/// @param dofs The dofs that are on this entity
/// @param dof_converter A DofIdConverter for the cell on which the entity
/// resides
/// @param is_dof_known An \c OverlayBasedIsDofKnownLookup that will be used
/// to query if a dof is unknown, ie needs to be interpolated.

template < class LAD, int DIM >
void LinearCellInterpolator< LAD, DIM >::interpolate_entity(
    int var, const Element< DataType, DIM > &element, const Entity &ent,
    const std::vector< int > &dofs, const DofIdConverter< LAD, DIM > &dof_converter,
    OverlayBasedIsDofKnownLookup< LAD, DIM > &is_dof_known) {
  const DofPartition< DataType, DIM > &dof_partition = element.get_space().dof();
  unsigned cell_dimension = element.get_cell().tdim();

  if (ent.tdim() < cell_dimension) {
    // check if a neighboring cell has already been
    // interpolated completely
    if (has_entity_interpolated_neighbor_cells(ent, element)) {
      // mark all dofs on this entity as known and return
      for (std::size_t i = 0; i < dofs.size(); ++i) {
        int local_id = -1;
        dof_converter.map_g2l(dofs[i], local_id);
        assert(local_id != -1);

        is_dof_known.mark_dof_as_known(local_id);
      }
      return;
    }
  }

  for (std::size_t i = 0; i < dofs.size(); ++i) {
    int local_id = -1;
    dof_converter.map_g2l(dofs[i], local_id);
    assert(local_id != -1);

    // overlay must be enabled to avoid interpolating the same
    // dof several times
    is_dof_known.enable_overlay();
    if (dof_partition.is_dof_on_subdom(dofs[i]) && !is_dof_known(dofs[i])) {
      is_dof_known.use_overlay(interpolate_from_interpolated_);

      unsigned num_contributions = 0;
      DataType sum_contributions = 0.0;

      std::vector< int > found_nearest_dofs;
      found_nearest_dofs.reserve(12);

      if (ent.tdim() < element.get_cell().tdim()) {

        // for all cells that contain dof
        for (IncidentEntityIterator cell = ent.begin_incident(cell_dimension);
             cell != ent.end_incident(cell_dimension); ++cell) {
          Element< DataType, DIM > neighboring_cell(element.get_space(), cell->index());

          DofIdConverter< LAD, DIM > neighboring_converter(&neighboring_cell);

          OverlayBasedIsDofKnownLookup< LAD, DIM > neighboring_lookup(
              dof_ident_, &neighboring_converter);

          neighboring_lookup.use_overlay(interpolate_from_interpolated_);

          OverlayBasedIsDofKnownLookup< LAD, DIM > *used_lookup =
              &neighboring_lookup;
          if (neighboring_cell.cell_index() == element.cell_index()) {
            // we are looking at the cell that we are currently interpolating,
            // we can just use our old is dof known lookup
            used_lookup = &is_dof_known;
          } else {
            // if we want to interpolate from interpolated dofs,
            // we need to adjust the is dof known lookup
            if (interpolate_from_interpolated_) {
              neighboring_lookup.resize_overlay_to_var(neighboring_cell, var);

              // configure IsDofKnownLookup such that all dofs
              // are known except the ones that are marked in
              // the still unknown bitmap
              for (int local_dof_id = 0;
                   local_dof_id < neighboring_cell.nb_dof(var);
                   ++local_dof_id) {
                if (!is_local_dof_still_unknown(
                        local_dof_id, var, neighboring_cell, dof_partition))
                  neighboring_lookup.mark_dof_as_known(local_dof_id);
              }
            }
          }

          NearestKnownDofsFinder< LAD, DIM > neighboring_finder(
              &neighboring_converter, used_lookup);

          std::vector< int > nearest_dofs;

          int neighbor_local_id = -1;
          if (neighboring_converter.map_g2l(dofs[i], neighbor_local_id)) {
            neighboring_finder.find_nearest_known_dofs(neighbor_local_id, var,
                                                       nearest_dofs);
            for (std::vector< int >::const_iterator it = nearest_dofs.begin();
                 it != nearest_dofs.end(); ++it) {
              // make sure dofs only contribute once
              if (std::find(found_nearest_dofs.begin(),
                            found_nearest_dofs.end(),
                            *it) == found_nearest_dofs.end()) {
                found_nearest_dofs.push_back(*it);
              }
            }
          }
        }
      } else // dim == cell tdim
      {
        NearestKnownDofsFinder< LAD, DIM > nearest_dof_finder(&dof_converter,
                                                         &is_dof_known);

        nearest_dof_finder.find_nearest_known_dofs(local_id, var,
                                                   found_nearest_dofs);
      }

      // std::vector<int> dofset_to_optimize = found_nearest_dofs;
      // if(found_nearest_dofs.size() != 12)
      //    subset_calculator(dofset_to_optimize, dofs[i], found_nearest_dofs);

      num_contributions = found_nearest_dofs.size();

      // 1 known Dof is not enough to interpolate..
      if (num_contributions > 1) {
        std::vector< DataType > nearest_values(num_contributions);

        vector_->GetValues(util::raw_array(found_nearest_dofs),
                           num_contributions, util::raw_array(nearest_values));

        sum_contributions =
            std::accumulate(nearest_values.begin(), nearest_values.end(), 0.0);

        DataType new_value = sum_contributions / num_contributions;

        vector_->SetValues(&dofs[i], 1, &new_value);

        is_dof_known.mark_dof_as_known(local_id);
        mark_dof_in_still_unknown_map(dof_partition, dofs[i], false);
      } else {
        mark_dof_in_still_unknown_map(dof_partition, dofs[i], true);
      }
    }
  }
}

/// Front end to interpolate a vector.
/// The interpolation only works if the mesh on which the interpolation takes
/// place is not the coarsest mesh in the hierarchy. The algorithm used here
/// can only interpolate vectors that have been transferred from the next
/// coarser grid.

template < class LAD, int DIM > class LinearInterpolation {
public:
  TYPE_FROM_CLASS(LAD, VectorType);

  LinearInterpolation< LAD, DIM >() : use_interpolation_matrix_(false) {}

  LinearInterpolation< LAD, DIM >(const bool use) : use_interpolation_matrix_(use) {}

  void use_interpolation_matrix(const bool use) {
    use_interpolation_matrix_ = use;
  }

  template < class ConnectedLevelType >
  void use_transposed_of_restriction_matrix(ConnectedLevelType &lvl) {
    if (lvl.is_scheduled_to_this_process()) {
      assert(lvl.restriction_matrix());

      lvl.create_interpolation_matrix();

      lvl.interpolation_matrix()->CreateTransposedFrom(
          *(lvl.restriction_matrix().get()));
    }
  }

  template < class ConnectedLevelType >
  void build_interpolation_matrix(ConnectedLevelType &lvl) {
    if (lvl.is_scheduled_to_this_process()) {
      std::cout << "\nERROR: LinearInterpolation::build_interpolation_matrix "
                   "not yet implemented.\n\n";
      exit(-1);
    }
  }

  /// Interpolates a vector
  /// @param level A \c ConnectedSingleLevel object in which the vector
  /// resides that shall be interpolated. The level supplied here must
  /// not be the coarsest level of the hierarchy.
  /// @param vector The vector that shall be interpolated

  template < class ConnectedLevelType >
  inline void operator()(const ConnectedLevelType &level, VectorType &vector) {
    if (level.is_scheduled_to_this_process()) {
      if (use_interpolation_matrix_) {
        assert(level.tmp_vector());
        assert(level.interpolation_matrix());

        level.interpolation_matrix()->VectorMult(vector,
                                                 level.tmp_vector().get());
        vector.CopyFrom(*(level.tmp_vector()));
      } else {
        vector.Update();

        unsigned dimension = level.mesh()->tdim();

        ForEachCell< LAD, DIM > for_each_cell(level.mesh(), level.space());

        assert(level.get_connection_to_next_coarser_grid() != NULL);

        // Bitmap to record which dofs could not be interpolated
        // after the first stage of the interpolation
        std::vector< bool > first_stage_unknown_dofs;

        std::vector< bool > *supplied_unknown_dofs_map = NULL;
        if (dimension == 3) {
          first_stage_unknown_dofs.resize(level.space()->dof().nb_dofs_global(),
                                          false);
          supplied_unknown_dofs_map = &first_stage_unknown_dofs;
        }

        LinearCellInterpolator< LAD, DIM > linear_cell_interpolator(
            level.get_connection_to_next_coarser_grid()
                ->get_dof_identification(),
            &vector, supplied_unknown_dofs_map, false);

        for_each_cell(linear_cell_interpolator);

        // Second interpolation to make sure all special cases in tetrahedrons
        // are covered
        if (dimension == 3) {
          LinearCellInterpolator< LAD, DIM > second_stage_linear_cell_interpolator(
              level.get_connection_to_next_coarser_grid()
                  ->get_dof_identification(),
              &vector, &first_stage_unknown_dofs, true);

          for_each_cell(second_stage_linear_cell_interpolator);
        }
      }
    }
  }

private:
  bool use_interpolation_matrix_;
};



} // namespace gmg
} // namespace la
} // namespace hiflow

#endif
