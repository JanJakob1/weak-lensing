set(space_SOURCES
  element.cc
  fe_evaluation.cc
  vector_space.cc
  )

# The following sources are needed for developing purpose only
if(EXPLICIT_TEMPLATE_INSTANTS)
  list(APPEND space_SOURCES   
    space.cc
  )
endif()

set(space_PUBLIC_HEADERS 
  dirichlet_boundary_conditions.h
  element.h
 # periodic_boundary_conditions.h
 # periodic_boundary_conditions_cartesian.h
  fe_evaluation.h
  fe_interpolation_cell.h
  fe_interpolation_global.h
  fe_interpolation_map.h
  vector_space.h
  )
