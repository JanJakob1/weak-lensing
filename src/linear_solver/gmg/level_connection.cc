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

#include "level_connection.h"


namespace hiflow {
namespace la {
namespace gmg {

template class BasicConnection< LADescriptorCoupledS,3 >;
template class SolutionTransferringConnection< LADescriptorCoupledS, 3 >;
template class GMGConnection< LADescriptorCoupledS, 3 >;
template class ConnectedLevel< LADescriptorCoupledS, 3, BasicConnection< LADescriptorCoupledS, 3 > >;
template class BasicHierarchy< ConnectedLevel< LADescriptorCoupledS, 3, BasicConnection< LADescriptorCoupledS, 3 > > >;

template class ConnectedLevel< LADescriptorCoupledS, 3, SolutionTransferringConnection< LADescriptorCoupledS, 3 > >;
    
template class GMGLevel< LADescriptorCoupledS, 3, GMGConnection< LADescriptorCoupledS, 3 > >;  
template class BasicHierarchy< GMGLevel< LADescriptorCoupledS, 3, GMGConnection< LADescriptorCoupledS, 3 > > >; 
      
template void visualize_multilevel_solutions<GMGLevel< LADescriptorCoupledD, 3, GMGConnection< LADescriptorCoupledD, 3 > > >(
    const BasicHierarchy< GMGLevel< LADescriptorCoupledD, 3, GMGConnection< LADescriptorCoupledD, 3 > > > &,
    const std::string &);

} // namespace gmg
} // namespace la
} // namespace hiflow
