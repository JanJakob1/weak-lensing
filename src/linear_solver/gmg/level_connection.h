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

#ifndef LEVEL_CONNECTION_H
#define LEVEL_CONNECTION_H

#include "vector_transfer.h"
#include "basic_hierarchy.h"
#include "gmg_level_impl.h"

namespace hiflow {
namespace la {
namespace gmg {

/// Manages all objects that exist between two levels in the Hierarchy:
/// Sets up DataTransferInformation, data transfer, DofIdentification.
/// The BasicConnection does not setup any transfer objects.

template < class LAD, int DIM > 
class BasicConnection {
public:
  /// Initializes the object. Collective on the communicator of the finer
  /// grid.
  /// @param coarse_level A coarser grid in the MultiLevelHierarchy
  /// @param fine_level A finer grid in the MultiLevelHierarchy.
  /// The communicator of the coarser grid must contain a subset of the
  /// processes in the communicator of the finer grid.
  BasicConnection(BasicLevel< LAD, DIM > &coarse_level,
                  BasicLevel< LAD, DIM > &fine_level);

  virtual ~BasicConnection() {}

  /// @return The DoF identification object

  boost::shared_ptr< const DofIdentification< LAD, DIM > >
  get_dof_identification() const {
    return dof_ident_;
  }

protected:
  boost::shared_ptr< DataTransferInformation< LAD, DIM > > info_;
  boost::shared_ptr< DofIdentification< LAD, DIM > > dof_ident_;

  boost::shared_ptr< DataTransferCoarseToFine< LAD, DIM > > transfer_to_fine_;
  boost::shared_ptr< DataTransferFineToCoarse< LAD, DIM > > transfer_to_coarse_;
};

/// An connection that allows the transfers of solution vectors

template < class LAD, int DIM >
class SolutionTransferringConnection : public BasicConnection< LAD, DIM > {
public:
  using BasicConnection< LAD, DIM >::dof_ident_;
  using BasicConnection< LAD, DIM >::info_;
  using BasicConnection< LAD, DIM >::transfer_to_coarse_;
  using BasicConnection< LAD, DIM >::transfer_to_fine_;

  /// Initializes the object. Collective on the communicator of the finer
  /// grid.
  /// @param coarse_level A coarser grid in the MultiLevelHierarchy
  /// @param fine_level A finer grid in the MultiLevelHierarchy.
  /// The communicator of the coarser grid must contain a subset of the
  /// processes in the communicator of the finer grid.

  SolutionTransferringConnection(BasicLevel< LAD, DIM > &coarse_level,
                                 BasicLevel< LAD, DIM > &fine_level);

  virtual ~SolutionTransferringConnection() {}

  /// @return A pointer to the vector transfer object that enables
  /// the transfer of the content between a vector on the coarse grid and
  /// a vector on the fine grid.

  boost::shared_ptr< VectorTransfer< LAD, DIM > > get_solution_transfer() const {
    return solution_transfer_;
  }

private:
  boost::shared_ptr< VectorTransfer< LAD, DIM > > solution_transfer_;
};

/// An connection that allows the transfers of solution vectors

template < class LAD, int DIM >
class GMGConnection : public SolutionTransferringConnection< LAD, DIM > {
public:
  typedef GMGLevel< LAD, DIM, GMGConnection< LAD, DIM > > LevelType;

  using BasicConnection< LAD, DIM >::dof_ident_;
  using BasicConnection< LAD, DIM >::info_;
  using BasicConnection< LAD, DIM >::transfer_to_coarse_;
  using BasicConnection< LAD, DIM >::transfer_to_fine_;

  /// Initializes the object. Collective on the communicator of the finer
  /// grid.
  /// @param coarse_level A coarser grid in the MultiLevelHierarchy
  /// @param fine_level A finer grid in the MultiLevelHierarchy.
  /// The communicator of the coarser grid must contain a subset of the
  /// processes in the communicator of the finer grid.

  GMGConnection(LevelType &coarse_level, LevelType &fine_level);

  virtual ~GMGConnection() {}

  /// Transfers the residual of the fine grid to the right hand side
  /// of the coarse grid

  void transfer_fine_res_to_coarse_rhs(void) {
    assert(fine_res_to_coarse_rhs_transfer_ != NULL);
    fine_res_to_coarse_rhs_transfer_->transfer_to_coarse();
  }

  /// Transfers the solution of the coarse grid to the residual of the fine
  /// grid.

  void transfer_coarse_sol_to_fine_res(void) {
    assert(coarse_sol_to_fine_res_transfer_ != NULL);
    coarse_sol_to_fine_res_transfer_->transfer_to_fine();
  }

  /// Transfers the solution of the coarse grid to the solution vector on
  /// the fine grid.

  void transfer_coarse_sol_to_fine_sol() {
    assert(this->get_solution_transfer() != NULL);
    this->get_solution_transfer()->transfer_to_fine();
  }

private:
  boost::shared_ptr< VectorTransfer< LAD, DIM > > fine_res_to_coarse_rhs_transfer_;
  boost::shared_ptr< VectorTransfer< LAD, DIM > > coarse_sol_to_fine_res_transfer_;
};

template < class LAD, int DIM >
BasicConnection< LAD, DIM >::BasicConnection(BasicLevel< LAD, DIM > &coarse_level,
                                        BasicLevel< LAD, DIM > &fine_level) {
  info_ = boost::shared_ptr< DataTransferInformation< LAD, DIM > >(
      new DataTransferInformation< LAD, DIM >(&coarse_level, &fine_level));

  transfer_to_fine_ = boost::shared_ptr< DataTransferCoarseToFine< LAD, DIM > >(
      new DataTransferCoarseToFine< LAD, DIM >(*info_));

  transfer_to_coarse_ = boost::shared_ptr< DataTransferFineToCoarse< LAD, DIM > >(
      new DataTransferFineToCoarse< LAD, DIM >(*info_));

  dof_ident_ = boost::shared_ptr< DofIdentification< LAD, DIM > >(
      new DofIdentification< LAD, DIM >(fine_level, coarse_level, info_,
                                   transfer_to_coarse_));

  dof_ident_->identify_dofs();
}

template < class LAD, int DIM >
SolutionTransferringConnection< LAD, DIM >::SolutionTransferringConnection(
    BasicLevel< LAD, DIM > &coarse_level, BasicLevel< LAD, DIM > &fine_level)
    : BasicConnection< LAD, DIM >(coarse_level, fine_level) {

  solution_transfer_ =
      boost::shared_ptr< VectorTransfer< LAD, DIM > >(new VectorTransfer< LAD, DIM >(
          dof_ident_, info_, fine_level.sol(), coarse_level.sol(),
          transfer_to_fine_, transfer_to_coarse_));
}

template < class LAD, int DIM >
GMGConnection< LAD, DIM >::GMGConnection(LevelType &coarse_level,
                                         LevelType &fine_level)
    : SolutionTransferringConnection< LAD, DIM >(coarse_level, fine_level) {

  fine_res_to_coarse_rhs_transfer_ =
      boost::shared_ptr< VectorTransfer< LAD, DIM > >(new VectorTransfer< LAD, DIM >(
          dof_ident_, info_, fine_level.res(), coarse_level.rhs(),
          transfer_to_fine_, transfer_to_coarse_));

  coarse_sol_to_fine_res_transfer_ =
      boost::shared_ptr< VectorTransfer< LAD, DIM > >(new VectorTransfer< LAD, DIM >(
          dof_ident_, info_, fine_level.res(), coarse_level.sol(),
          transfer_to_fine_, transfer_to_coarse_));
}

} // namespace gmg
} // namespace la
} // namespace hiflow

#endif
