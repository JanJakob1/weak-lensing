set(gmg_SOURCES
  # these are need for compiling tests only
  gmg_level.cc
  data_transfer.cc
  #dof_identification.cc
  vector_transfer.cc
  #level_connection.cc
  #basic_hierarchy.cc		#TODO> LocalObjectType local_obj’ has incomplete type typename AssemblyPolicy::LocalObjectType local_obj;
  gmg_hierarchy.cc
  #gmg_interpolation.cc
  #gmg_restriction.cc
  gmg.cc
  auxiliary.cc
)

set(gmg_PUBLIC_HEADERS
  util.h		#should work
  hybrid_map.h
  gmg_level.h
  data_transfer.h
  dof_identification.h
  vector_transfer.h
  level_connection.h
  basic_hierarchy.h
  gmg_hierarchy.h
  gmg_interpolation.h
  gmg_restriction.h
  gmg.h
  gmg_level_impl.h
)
