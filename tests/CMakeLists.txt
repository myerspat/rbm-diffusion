# Add new header files here
set(HEADER_FILES
  unit_test_framework.hpp
)

# Add new source files here
set(SOURCE_FILES
  test_matrix.cpp
  test_linalg.cpp
  test_mesh.cpp
)

set(LIBRARY_NAME cpplib)

foreach(test ${SOURCE_FILES})
  set(test_name ${test}_ctest)
  add_executable(${test_name} ${test})
  target_link_libraries(${test_name} cpplib)
  add_test(${test_name} ${test_name})
endforeach()
