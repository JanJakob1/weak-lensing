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

#include "gmg_level.h"


namespace hiflow {
namespace la {
namespace gmg {

template class BasicLevel< LADescriptorCoupledS, 3 >;
template < class ConnectionType > class ConnectedLevel< LADescriptorCoupledS, 3, ConnectionType >;
template < class ConnectionType > class GMGLevel< LADescriptorCoupledS, 3, ConnectionType >;

} // namespace gmg
} // namespace la
} // namespace hiflow
