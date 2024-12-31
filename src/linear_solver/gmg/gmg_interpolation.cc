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

#include "gmg_interpolation.h"

namespace hiflow {
namespace la {
namespace gmg {




/// For a given cell and variable, returns an ijk iterator of corresponding type


template class IjkIteratorFactory< LADescriptorCoupledD,3 >;
template class IjkIteratorFactory< LADescriptorCoupledS, 3 >;

template class DofIdConverter< LADescriptorCoupledD,3 >;
template class DofIdConverter< LADescriptorCoupledS,3 >;

template class OnDemandDofIdConverter< LADescriptorCoupledD,3 >;
template class OnDemandDofIdConverter< LADescriptorCoupledS,3 >;

template class IsDofKnownLookup< LADescriptorCoupledD,3 >;
template class IsDofKnownLookup< LADescriptorCoupledS ,3>;

template class DirectIsDofKnownLookup< LADescriptorCoupledD ,3>;
template class DirectIsDofKnownLookup< LADescriptorCoupledS,3 >;

template class OverlayBasedIsDofKnownLookup< LADescriptorCoupledD,3 >;
template class OverlayBasedIsDofKnownLookup< LADescriptorCoupledS,3 >;



/// Finds the nearest known dof on a subentity
/// @return The global dof id of the nearest known dof or -1, if no
/// known dof has been found.
/// @param tdim The dimension of the subentity
/// @param sindex The sindex of the subentity (unused if tdim == tdim of cell)
/// @param global_dof The global dof of which the nearest dofs shall be found.
/// It must also reside on the supplied subentity.

/// Performs the search for nearest known dofs on Q (quad based) cells
/// @param origin A ijk dof id representing the origin of the search
/// @param var The variable to which the ijk dof id \c origin belongs
/// @param found_dofs A vector that will contain the nearest known dofs
/// that have been found once a call to this function returns.

template class NearestKnownDofsFinder< LADescriptorCoupledD,3 >;
template class NearestKnownDofsFinder< LADescriptorCoupledS,3 >;

template class ForEachCell< LADescriptorCoupledD,3 >;
template class ForEachCell< LADescriptorCoupledS,3 >;

/// Interpolates the unkown dofs on a cell.
/// This only works if the mesh on which the interpolation takes place is
/// not the coarsest mesh in the hierarchy. The algorithm used here can only
/// interpolate vectors that have been transferred from the next coarser grid.


/// Interpolate an entity
/// @param var The variable id
/// @param element The cell on which the entity resides
/// @param ent The entity that shall be interpolated
/// @param dofs The dofs that are on this entity
/// @param dof_converter A DofIdConverter for the cell on which the entity
/// resides
/// @param is_dof_known An \c OverlayBasedIsDofKnownLookup that will be used
/// to query if a dof is unknown, ie needs to be interpolated.

template class LinearCellInterpolator< LADescriptorCoupledD,3 >;
template class LinearCellInterpolator< LADescriptorCoupledS,3 >;

template class LinearInterpolation< LADescriptorCoupledD,3 >;
template class LinearInterpolation< LADescriptorCoupledS,3 >;

} // namespace gmg
} // namespace la
} // namespace hiflow
