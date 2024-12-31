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

#include "dof_identification.h"
#include "gmg_level.h"
#include "mesh/iterator.h"
#include "mesh/mesh.h"

namespace hiflow {
namespace la {
namespace gmg {

void do_nothing(int unused) {}

namespace detail {
/// Compares two exchanged_values objects,
/// and enables the use of exchanged_values objects as key in
/// boost::unordered_map hashtables. The rank stored in the exchanged_values
/// objects is not compared, because two exchanged_values objects shall be equal
/// if their coordinates and variable ids are equal. On which process the DoFs
/// are located is not relevant, because we want to know which DoFs correspond
/// to each other.
/// @return whether the DoFs tow which the exchanged_values objects refer
/// correspond to each other, i.e. whether the coordinates and variable ids
/// of the exchanged_values objects are equal.
/// @param lhs left hand side exchanged_values object
/// @param rhs right hand side exchanged_values object

bool operator==(const DofIdentificationExchangedValues &lhs,
                const DofIdentificationExchangedValues &rhs) {
  double rel_epsilon = 100.0 * std::numeric_limits< double >::epsilon();

  if (lhs.coordinates.size() != rhs.coordinates.size())
    return false;

  for (std::size_t i = 0; i < lhs.coordinates.size(); ++i)
    if (std::abs(lhs.coordinates[i] - rhs.coordinates[i]) >
        rel_epsilon * std::abs(lhs.coordinates[i]))
      return false;

  if (lhs.variable != rhs.variable)
    return false;
  // We do not check for the rank, because where it lies is not
  // important for the equality of two dofs

  return true;
}

/// Hashes an exchanged_values object.
/// This enables the use of exchanged_values objects as key in
/// boost::unordered_map hashtables, which optimizes the coordinate comparison
/// of two DoF lists. (see
/// DofIdentification::generate_coordinate_intersection_map(...)) The rank
/// stored in an exchanged_values object is not hashed, because two
/// exchanged_values objects shall be equal if their coordinates and variable
/// ids are equal. On which process the DoFs are located is not relevant,
/// because we want to know which DoFs correspond to each other.
/// @return The hashed value for an exchanged_values object.
/// @param data The exchanged_values object to hash.
std::size_t hash_value(const DofIdentificationExchangedValues &data) {
  boost::hash< std::vector< float > > coordinate_hasher;
  std::vector< float > float_coordinates(data.coordinates.size());
  for (std::size_t i = 0; i < data.coordinates.size(); ++i)
    float_coordinates[i] = static_cast< float >(data.coordinates[i]);

  std::size_t result = coordinate_hasher(float_coordinates);

  boost::hash< int > int_hasher;

  // We do not hash the rank
  result ^= int_hasher(data.variable);

  return result;
}

}
} // namespace gmg
} // namespace la
} // namespace hiflow
