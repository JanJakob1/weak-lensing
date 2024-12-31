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

/// \author Aksel Alpay, Martin Wlotzka

#include "gmg_hierarchy.h"
#include "basic_hierarchy.h"
#include "data_transfer.h"
#include "dof_identification.h"
#include "gmg_level.h"
#include "vector_transfer.h"
#include "level_connection.h"
#include "gmg_restriction.h"
#include "gmg_interpolation.h"
#include "gmg.h"

namespace hiflow {
namespace la {
namespace gmg {

template class BasicHierarchy< BasicLevel< LADescriptorCoupledS, 3 > >;
template void visualize_multilevel_solutions< BasicLevel< LADescriptorCoupledS, 3 > >(
                const BasicHierarchy< BasicLevel< LADescriptorCoupledS, 3 > > &,
                const std::string &);

template class HierarchyTypes< LADescriptorCoupledS, 3 >;

template class DataTransferFineToCoarse< LADescriptorCoupledS, 3 >;
template class DataTransferCoarseToFine< LADescriptorCoupledS, 3 >;
template class DataTransferInformation< LADescriptorCoupledS, 3 >;

template class DofIdentification< LADescriptorCoupledS, 3 >;

template class BasicLevel< LADescriptorCoupledS, 3 >;
template < class ConnectionType > class ConnectedLevel< LADescriptorCoupledS, 3, ConnectionType >;
template < class ConnectionType > class GMGLevel< LADescriptorCoupledS, 3, ConnectionType >;

template class VectorTransfer< LADescriptorCoupledS, 3 >;

template class BasicConnection< LADescriptorCoupledS,3 >;
template class SolutionTransferringConnection< LADescriptorCoupledS, 3 >;
//template class GMGConnection< LADescriptorCoupledS, 3 >;
//template class ConnectedLevel< LADescriptorCoupledS, 3, BasicConnection< LADescriptorCoupledS, 3 > >;
template class BasicHierarchy< ConnectedLevel< LADescriptorCoupledS, 3, BasicConnection< LADescriptorCoupledS, 3 > > >;

template class ConnectedLevel< LADescriptorCoupledS, 3, SolutionTransferringConnection< LADescriptorCoupledS, 3 > >;
    
//template class GMGLevel< LADescriptorCoupledS, 3, GMGConnection< LADescriptorCoupledS, 3 > >;  
template class BasicHierarchy< GMGLevel< LADescriptorCoupledS, 3, GMGConnection< LADescriptorCoupledS, 3 > > >; 
      
template void visualize_multilevel_solutions<GMGLevel< LADescriptorCoupledD, 3, GMGConnection< LADescriptorCoupledD, 3 > > >(
    const BasicHierarchy< GMGLevel< LADescriptorCoupledD, 3, GMGConnection< LADescriptorCoupledD, 3 > > > &,
    const std::string &);

template class LinearCellRestrictor< LADescriptorCoupledS, 3 >;
template class LinearRestriction< LADescriptorCoupledS, 3 >;

template class DofIdConverter< LADescriptorCoupledS, 3 >;
template class OnDemandDofIdConverter< LADescriptorCoupledS, 3 >;
template class IsDofKnownLookup< LADescriptorCoupledS, 3 >;
template class DirectIsDofKnownLookup< LADescriptorCoupledS, 3 >;
template class OverlayBasedIsDofKnownLookup< LADescriptorCoupledS, 3 >;
template class NearestKnownDofsFinder< LADescriptorCoupledS, 3 >;
template class ForEachCell< LADescriptorCoupledS, 3 >;
template class IjkIteratorFactory< LADescriptorCoupledS, 3 >;
template class LinearCellInterpolator< LADescriptorCoupledS, 3 >;
template class LinearInterpolation< LADescriptorCoupledS, 3 >;

template class SmootherHierarchy< LADescriptorCoupledD >;
template class SmootherHierarchy< LADescriptorCoupledS >;

template class SetupJacobiSmoother< LADescriptorCoupledD >;
template class SetupJacobiSmoother< LADescriptorCoupledS >;

template class SetupAsyncIterGPUSmoother< LADescriptorCoupledD >;
template class SetupAsyncIterGPUSmoother< LADescriptorCoupledS >;

template class SetupIndividualSmoother< LADescriptorCoupledD >;
template class SetupIndividualSmoother< LADescriptorCoupledS >;

} // namespace gmg

template class GeometricMultiGrid< LADescriptorCoupledS, 3 >;

} // namespace la
} // namespace hiflow
