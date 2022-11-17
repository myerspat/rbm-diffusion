# Include path search paths
set(OPENBLAS_INCLUDE_SEARCH_PATHS
  /usr/include
  /usr/include/openblas
  /usr/local/include
  /usr/local/include/openblas
  /opt/OpenBLAS/include
  $ENV{OpenBLAS_HOME}
  $ENV{OpenBLAS_HOME}/include
)

# Library search paths
set(OPENBLAS_LIB_SEARCH_PATHS
  /lib
  /lib/openblas
  /usr/lib
  /usr/local/lib
  /opt/OpenBLAS/lib
  $ENV{OpenBLAS}cd
  $ENV{OpenBLAS}/lib
  $ENV{OpenBLAS_HOME}
  $ENV{OpenBLAS_HOME}/lib
)

# Find path and library
find_path(OPENBLAS_INCLUDE_DIR NAMES cblas.h PATHS ${OPENBLAS_INCLUDE_SEARCH_PATHS})
find_library(OPENBLAS_LIB NAMES openblas PATHS ${OPENBLAS_LIB_SEARCH_PATHS})

# Print status
set(OpenBLAS_FOUND ON)
if(NOT OPENBLAS_INCLUDE_DIR)
  set(OpenBLAS_FOUND OFF)
  message(STATUS "Could not find OpenBLAS include")
endif()
if(NOT OPENBLAS_LIB)
  set(OpenBLAS_FOUND OFF)
  message("Could not find OpenBLAS library")
endif()

# Output paths if found, throw errors if required but not found
if(OpenBLAS_FOUND)
  if(NOT OpenBLAS_FIND_QUIETLY)
    message(STATUS "Found OpenBLAS libraries: ${OPENBLAS_LIB}")
    message(STATUS "Found OpenBLAS include: ${OPENBLAS_INCLUDE_DIR}")
  endif()
else()
  if(OpenBLAS_FIND_REQUIRED)
    message(FATAL_ERROR "Could not find OpenBLAS")
  endif()
endif()

mark_as_advanced(
  OPENBLAS_INCLUDE_DIR
  OPENBLAS_LIB
  OpenBLAS
)
