set(cell_trafo_SOURCES
)

# the following source files are only included for developping purpose
if(EXPLICIT_TEMPLATE_INSTANTS)
  list(APPEND cell_trafo_SOURCES   
    cell_trafo.cc
)
endif()

set(cell_trafo_PUBLIC_HEADERS
 cell_transformation.h
 linearlinetransformation.h
 lineartriangletransformation.h
 lineartetrahedrontransformation.h
 linearpyramidtransformation.h
 bilinearquadtransformation.h
 trilinearhexahedrontransformation.h
 alignedquadtransformation.h
 alignedhexahedrontransformation.h
 )
