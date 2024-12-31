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

#include "nonlinear/damping_armijo.h"

namespace hiflow {

/// template instantiation
<<<<<<< HEAD
template class ArmijoDamping< la::LADescriptorCoupledD, 3 >;

/*
template class ArmijoDamping< la::LADescriptorCoupledS, 3 >;
#ifdef WITH_HYPRE
template class ArmijoDamping< la::LADescriptorHypreD, 3 >;
#endif
template class ArmijoDamping< la::LADescriptorPolynomialChaosD, 3 >;
#ifdef WITH_HYPRE
template class ArmijoDamping< la::LADescriptorPolynomialChaosExpansionD, 3 >;
#endif
=======
template class ArmijoDamping< la::LADescriptorCoupledD >;
template class ArmijoDamping< la::LADescriptorCoupledS >;
template class ArmijoDamping<
    la::LADescriptorPCE< la::LADescriptorCoupledD > >;
template class ArmijoDamping<
    la::LADescriptorPCE< la::LADescriptorCoupledS > >;
>>>>>>> master
template class ArmijoDamping<
    la::LADescriptorBlock< la::LADescriptorCoupledD >, 3 >;
template class ArmijoDamping<
<<<<<<< HEAD
    la::LADescriptorBlock< la::LADescriptorCoupledS >, 3 >;
#ifdef WITH_HYPRE
template class ArmijoDamping< la::LADescriptorBlock< la::LADescriptorHypreD >, 3 >;
=======
    la::LADescriptorBlock< la::LADescriptorCoupledS > >;
template class ArmijoDamping< la::LADescriptorPolynomialChaosD >;

#ifdef WITH_HYPRE
template class ArmijoDamping< la::LADescriptorHypreD >;
template class ArmijoDamping< la::LADescriptorPCE< la::LADescriptorHypreD > >;
template class ArmijoDamping< la::LADescriptorBlock< la::LADescriptorHypreD > >;
>>>>>>> master
#endif
*/
} // namespace hiflow
