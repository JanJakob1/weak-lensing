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

#include "degree_of_freedom.h"

/// \author Michael Schick<br>Martin Baumann

namespace hiflow {
namespace doffem {

template class DegreeOfFreedom< float, 3 >;
template void permute_constrained_dofs_to_end(DegreeOfFreedom< float, 3 > &dof);

} // namespace doffem
} // namespace hiflow
