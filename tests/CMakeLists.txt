# Add new header files here
set(HEADER_FILES
  unit_test_framework.hpp
)

# Add new source files here
set(TEST_FILES
  test_mesh
  test_input
  test_meshElement
  test_cell
  test_rbm
  test_material
)

foreach(test ${TEST_FILES})
  add_executable(${test} ${test}.cpp)
  target_link_libraries(${test} rbmlib)
  add_test(${test} ${test})
endforeach()
