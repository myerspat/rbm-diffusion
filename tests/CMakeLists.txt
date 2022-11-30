# Add new header files here
set(HEADER_FILES
  unit_test_framework.hpp
)

# Add new source files here
set(TEST_FILES
  test_mesh
  test_input
  test_meshElement
  test_rbm
  test_material
)

# For each test in TEST_FILES
foreach(test ${TEST_FILES})
  # Add test executable
  add_executable(${test} ${test}.cpp)

  # Link the test to rbm library
  target_link_libraries(${test} rbm)

  # Add test to ctest suite
  add_test(${test} ${test})
endforeach()

# Regression tests
set(TEST_FILES
  test_1D_mesh
)

# For each test in TEST_FILES
foreach(test ${TEST_FILES})
  # Add test executable
  add_executable(regression_${test} regression_tests/${test}.cpp)

  # Link the test to rbm library
  target_link_libraries(regression_${test} rbm)

  # Add test to ctest suite
  add_test(regression_${test} regression_${test})
endforeach()
