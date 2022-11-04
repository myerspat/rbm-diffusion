# Add new header files here
set(HEADER_FILES
  matrix.hpp
  linalg.hpp
  mesh.hpp
)

set(SOURCE_FILES
  matrix.cpp
  linalg.cpp
  mesh.cpp
)

add_library(cpplib ${SOURCE_FILES} ${HEADER_FILES})
