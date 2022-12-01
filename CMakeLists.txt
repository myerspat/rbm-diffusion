# Require CMAKE 3.1 of higher
cmake_minimum_required(VERSION 3.10.0)

# Project name
project(rbm VERSION 0.1.0)

# C++ supported version minimum
enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#================================================================================
# Compiler options

# Supported options
option(RBM_ENABLE_PROFILE "Compile with profiling flags"  OFF)
option(RBM_OPTIMIZED      "Compile with optimization"     OFF)
option(RBM_GIT_SUBMODULE  "Check submodules during build" ON)

# Profiling flags
if(RBM_ENABLE_PROFILE)
  list(APPEND flags -g -fno-omit-frame-pointer)
endif()

# Optimization flags
if(RBM_OPTIMIZED)
  list(APPEND flags -O3 -march=native)
endif()

# Show flags being used
message(STATUS "RBM C++ flags: ${flags}")

#================================================================================
# RBM library

# Source files
list(APPEND rbm_SOURCES
  src/mesh.cpp
  src/input.cpp
  src/material.cpp
  src/meshElement.cpp
  src/rbm.cpp
)

# Archived source files
list(APPEND rbm_SOURCES
  src/archive/matrix.cpp
  src/archive/linalg.cpp
)

# Include headers
include_directories(include)

# Add src to library
add_library(rbm SHARED ${rbm_SOURCES})

# Add compiler options
target_compile_options(rbm PRIVATE ${flags})

#================================================================================
# Submodules

# Check git and add git submodules if specified
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  # Update submodules if needed
  if(RBM_GIT_SUBMODULE)
    message(STATUS "Submodule update")
    execute_process(COMMAND git submodule update --init -- vendor/
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
    endif()
  endif()
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/vendor/xtensor/CMakeLists.txt")
  message(FATAL_ERROR "The submodules were not downloaded! RBM_GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()

# xtl library
if(NOT xtl)
  add_subdirectory(vendor/xtl)
  set(xtl_DIR ${CMAKE_CURRENT_BINARY_DIR}/vendor/xtl)
  message(STATUS "xtl module: ${xtl_DIR}")
endif()

# xtensor library
if(NOT xtensor)
  add_subdirectory(vendor/xtensor)
  set(xtensor_DIR ${CMAKE_CURRENT_BINARY_DIR}/vendor/xtensor)
  message(STATUS "xtensor module: ${xtensor_DIR}")
endif()

# xtensor-blas library
if(NOT xtensor-blas)
  add_subdirectory(vendor/xtensor-blas)
  set(xtensor-blas_DIR ${CMAKE_CURRENT_BINARY_DIR}/vendor/xtensor-blas)
  message(STATUS "xtensor-blas module: ${xtensor-blas_DIR}")
endif()

# pugixml library
if(NOT pugixml)
  add_subdirectory(vendor/pugixml)
  set(pugixml_DIR ${CMAKE_CURRENT_BINARY_DIR}/vendor/pugixml)
  message(STATUS "pugixml module: ${pugixml_DIR}")
endif()

#================================================================================
# Dependencies

# Include OpenBLAS library
include(cmake/FindOpenBLAS.cmake)

#================================================================================
# Linking / changing project executable name

set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "rbm-diffusion")
target_include_directories(rbm PRIVATE ${CMAKE_BINARY_DIR}/include)
target_link_libraries(rbm xtensor xtensor-blas pugixml ${OPENBLAS_LIB})

# Main executable
add_executable(rbm-diffusion src/main.cpp)
target_link_libraries(rbm-diffusion rbm)

#================================================================================
# CTest and subdirectories

include(CTest)
add_subdirectory(tests)
add_subdirectory(examples)
enable_testing()
