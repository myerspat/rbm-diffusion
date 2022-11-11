# Require CMAKE 3.1 of higher
cmake_minimum_required(VERSION 3.10.0)

# Project name
project(rbm VERSION 0.1.0)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(RBM_ENABLE_PROFILE "Compile with profiling flags" OFF)
option(RBM_OPTIMIZED      "Compile with optimization"    OFF)

if(RBM_ENABLE_PROFILE)
  list(APPEND flags -g -fno-omit-frame-pointer)
endif()

if(RBM_OPTIMIZED)
  list(APPEND flags -O3 -march=native)
endif()

# Show flags being used
message(STATUS "RBM C++ flags: ${flags}")

# Source files
list(APPEND rbmlib_SOURCES
  src/mesh.cpp
)

# Archived source files
list(APPEND rbmlib_SOURCES
  src/archive/matrix.cpp
  src/archive/linalg.cpp
)

include_directories(include)

add_library(rbmlib SHARED ${rbmlib_SOURCES})

target_compile_options(rbmlib PRIVATE ${flags})

# xtensor header-only library
add_subdirectory(vendor/xtl)
set(xtl_DIR ${CMAKE_CURRENT_BINARY_DIR}/vendor/xtl)
add_subdirectory(vendor/xtensor)
add_subdirectory(vendor/xtensor-blas)

# OpenBLAS
# set(BLA_VENDER OpenBLAS)
# find_package(BLAS REQUIRED)
# if(BLAS_FOUND)
#   message("-- OpenBLAS found: ${BLAS_LIBRARIES}")
#   target_include_directories(rbmlib PRIVATE ${CMAKE_BINARY_DIR}/opt/OpenBLAS/include/)
# endif(BLAS_FOUND)

target_include_directories(rbmlib PRIVATE ${CMAKE_BINARY_DIR}/opt/OpenBLAS/include/)

target_include_directories(rbmlib PRIVATE ${CMAKE_BINARY_DIR}/include)
target_link_libraries(rbmlib xtensor xtensor-blas /opt/OpenBLAS/lib/libopenblas.so)

include(CTest)
add_subdirectory(tests)
add_subdirectory(examples)

enable_testing()
