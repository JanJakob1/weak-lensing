project(HIFLOW CXX)

cmake_minimum_required(VERSION 3.0.2)

# Prevent in-source builds
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR "In-source builds are not allowed. Please create a build directory and run cmake from there, passing the path to this source directory as the last argument.
Before building, delete the file CMakeCache.txt and the directory CMakeFiles that were created in the source directory.")
endif()

# guard against bad build-type strings
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
message("Build type: ${CMAKE_BUILD_TYPE}")

# add cmake directory to module path to use our own modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

include(Version)

# Verbose compiler printing
set(CMAKE_VERBOSE_MAKEFILE false)

# Give hints about where to search for includes and libraries.
set(CMAKE_INCLUDE_PATH "/opt/intel/mkl/*/include;/opt/intel/mkl/10.2.5.035/include;/usr/include;/usr/include/mkl")
set(CMAKE_LIBRARY_PATH "/opt/intel/mkl/*/lib/em64t;/opt/intel/mkl/10.2.5.035/lib/em64t;/usr/lib;/usr/lib/x86_64-linux-gnu")

# Configure required and optional libraries, and generate config.h file.
include(HiFlowConfiguration)

# Set internal variables
set(HIFLOW_LIB_DIR lib)
set(HIFLOW_BIN_DIR bin)
set(HIFLOW_INCLUDE_DIR include/hiflow3)
set(HIFLOW_EXAMPLES_DIR share/hiflow3/examples)

# DEBUG / RELEASE build flags
add_definitions(-DMPICH_IGNORE_CXX_SEEK)  # Some MPI implementations need this.

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

set(CMAKE_CXX_FLAGS_COVERAGE "-coverage" CACHE STRING
    "Flags used by the C++ compiler during coverage builds."
    FORCE)
set(CMAKE_C_FLAGS_COVERAGE "-coverage" CACHE STRING
    "Flags used by the C compiler during coverage builds."
    FORCE)
set(CMAKE_EXE_LINKER_FLAGS_COVERAGE "-coverage" CACHE STRING
    "Flags used for linking binaries during coverage builds."
    FORCE)
set(CMAKE_SHARED_LINKER_FLAGS_COVERAGE "-coverage" CACHE STRING
    "Flags used by the shared libraries linker during coverage builds."
    FORCE)
mark_as_advanced(
    CMAKE_CXX_FLAGS_COVERAGE
    CMAKE_C_FLAGS_COVERAGE
    CMAKE_EXE_LINKER_FLAGS_COVERAGE
    CMAKE_SHARED_LINKER_FLAGS_COVERAGE )
# Update the documentation string of CMAKE_BUILD_TYPE for GUIs
set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel Coverage."
    FORCE)

# VENDOR-specific compiler flags
if (CMAKE_CXX_COMPILER MATCHES g[+][+].* OR CMAKE_CXX_COMPILER MATCHES c[+][+].*)
  message("Setting vendor-specific compiler flags for GNU Compiler")
  set(HIFLOW_VENDOR_FLAGS "-pthread -ansi")
endif()

if (CMAKE_CXX_COMPILER MATCHES icpc.*)
  message("Setting vendor-specific compiler flags for Intel Compiler")
  # Use interprocedural optimizations -- requires xiar archiver.
  set(HIFLOW_VENDOR_FLAGS "-ipo-separate")
  set(CMAKE_AR "xiar")
endif()

if (CMAKE_CXX_COMPILER MATCHES pgCC.*)
  message("Setting vendor-specific compiler flags for PGI Compiler")
  set(HIFLOW_VENDOR_FLAGS "-noswitcherror -cppsuffix=.cu")
endif()

if (CMAKE_CXX_COMPILER MATCHES clang[+][+].*)
  message("Setting vendor-specific compiler flags for Clang Compiler")
  set(HIFLOW_VENDOR_FLAGS "")
endif()

# Option for c++ iso standard
set(HIFLOW_CPP_ISO_STANDARD "c++11" CACHE STRING "C++ ISO standard")
if (HIFLOW_CPP_ISO_STANDARD STREQUAL "c++17")
  message("Use c++17 iso standard")
  set(CMAKE_CXX_STANDARD 17)
elseif (HIFLOW_CPP_ISO_STANDARD STREQUAL "c++14")
  message("Use c++14 iso standard")
  set(CMAKE_CXX_STANDARD 14)
elseif (HIFLOW_CPP_ISO_STANDARD STREQUAL "c++11")
  message("Use c++11 iso standard")
  set(CMAKE_CXX_STANDARD 11)
endif()

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${HIFLOW_VENDOR_FLAGS} -Wall -Wextra -Wno-reorder -Wno-unknown-pragmas -Wno-sign-compare -Wno-cast-function-type -Wno-comment -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable ")

# Boost include directories
include_directories(${Boost_INCLUDE_DIR})

# TinyXML include directory
include_directories(${TINYXML2_INCLUDE_DIR})

# Process subdirectories
add_subdirectory(src)
add_subdirectory(doc)
add_subdirectory(contrib)

if (BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

if (BUILD_UTILS)
  add_subdirectory(utils)
endif()

if (BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if (BUILD_EXERCISES)
  add_subdirectory(exercises)
endif()

# Install export file
install(EXPORT HiFlow DESTINATION share/hiflow3)

# Set these variable _before_ including CPack (why is this necessary???)
set(CPACK_INSTALL_CMAKE_PROJECTS ${PROJECT_BINARY_DIR} hiflow ALL /)

# CPack configuration variables
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "HiFlow3")
set(CPACK_PACKAGE_FILE_NAME "HiFlow3")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "3")
set(CPACK_PACKAGE_VERSION_PATCH "2")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${PROJECT_SOURCE_DIR}/doc/INSTALLATION")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/COPYING.GPL")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/COPYING.LGPL")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/AUTHORS")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "HiFlow3_version ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}.${CMake_VERSION_PATCH}")

include(CPack)
