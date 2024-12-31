// Copyright (C) 2011-2020 Vincent Heuveline
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

#include "gmg.h"

#include <csignal>

namespace hiflow {
namespace la {
namespace gmg {

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
