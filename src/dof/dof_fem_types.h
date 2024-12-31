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

#ifndef __DOF_FEM_TYPES_H_
#define __DOF_FEM_TYPES_H_
#include <memory>

/// \author Michael Schick<br>Martin Baumann

namespace hiflow {
namespace doffem {
  template <class DataType, int DIM> class CellTransformation;
  template <class DataType, int DIM> class RefCell;
  
  
  template <class DataType, int DIM> 
  using CellTrafoPtr = std::shared_ptr< CellTransformation<DataType, DIM> >;
  
  template <class DataType, int DIM> 
  using ConstCellTrafoPtr = std::shared_ptr< const CellTransformation<DataType, DIM> >;
  
  template <class DataType, int DIM> 
  using RefCellPtr = std::shared_ptr< RefCell<DataType, DIM> >;

  template <class DataType, int DIM> 
  using ConstRefCellPtr = std::shared_ptr< const RefCell<DataType, DIM> >;  
  
  typedef int DofID;
  typedef int FETypeID;
  
  enum RefCellType
  {
    REF_CELL_NOT_SET = 0,
    REF_CELL_LINE_STD = 1,
    REF_CELL_TRI_STD = 2,
    REF_CELL_QUAD_STD = 3,
    REF_CELL_TET_STD = 4,
    REF_CELL_HEX_STD = 5,
    REF_CELL_PYR_STD = 6,
  };

  // Type of dof collection defining a FE on reference cell
  enum DofContainerType
  {
    DOF_CONTAINER_NOT_SET = 0,
    DOF_CONTAINER_LAGRANGE = 1,
    DOF_CONTAINER_BDM = 2,
    DOF_CONTAINER_RT = 3
  };
  
  enum DofFunctionalType 
  {
    DOF_NOT_SET = 0,
    DOF_POINT_EVAL = 1,
    DOF_POINT_NORMAL_EVAL = 2,
    DOF_POINT_TANGENT_EVAL = 3,
    DOF_CELL_MOMENT = 4,
    DOF_FACET_MOMENT = 5
  };
  
  enum DofConstraint
  {
    DOF_CONSTRAINT_INTERIOR = 0,
    DOF_CONSTRAINT_FACET = 1,
    DOF_CONSTRAINT_EDGE = 2,
    DOF_CONSTRAINT_VERTEX = 3  
  };

  enum DofPosition
  {
    DOF_POSITION_INTERIOR = 0,
    DOF_POSITION_FACET = 1,
    DOF_POSITION_EDGE = 2,
    DOF_POSITION_VERTEX = 3,  
  };
  
  /// Enumeration of different DoF ordering strategies. HIFLOW_CLASSIC refers to
  /// the DoF numbering as always done in HiFlow3. The other two options allow
  /// permutations of the classic numbering by means of the Cuthill-McKee and the
  /// King method, respectively.

  enum DOF_ORDERING 
  { 
    HIFLOW_CLASSIC, 
    CUTHILL_MCKEE, 
    KING 
  };

  enum AnsatzSpaceType 
  {
    ANSATZ_NOT_SET = 0,
    ANSATZ_P_LAGRANGE = 1,
    ANSATZ_Q_LAGRANGE = 2,
    ANSATZ_SKEW_P_AUG = 3,
    ANSATZ_P_AUG = 4,
    ANSATZ_RT = 5,
    ANSATZ_SUM = 6
  };
  
  // FE type on reference cell
  enum FEType 
  {
    FE_TYPE_NOT_SET = 0,
    FE_TYPE_CUSTOM = 1,
    FE_TYPE_LAGRANGE = 2,
    FE_TYPE_BDM = 3,
    FE_TYPE_RT = 4,
    FE_TYPE_LAGRANGE_VECTOR = 5
  };
 
  enum FEConformity
  {
    FE_CONFORM_NONE = -1,
    FE_CONFORM_L2 = 0,
    FE_CONFORM_HDIV = 1,
    FE_CONFORM_HCURL = 2,
    FE_CONFORM_H1 = 3
  };

  enum FETransformationType
  {
    FE_TRAFO_STD = 0,
    FE_TRAFO_PIOLA = 1
  };
} // namespace doffem
} // namespace hiflow

#endif
