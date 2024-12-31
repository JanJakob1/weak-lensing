# get the main sources
file(GLOB SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} lib/*.cpp lib/*.h)
source_group("" FILES ${SRCS})

set(Version 1.1.1)

# create the lib
add_library(iluplusplus-${Version} STATIC ${SRCS})
target_compile_options(iluplusplus-${Version} PUBLIC -Wno-write-strings -Wno-format-security -Wno-unused-result)
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/lib/" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/lib/")
