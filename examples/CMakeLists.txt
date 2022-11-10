# Set example files
set(EXAMPLES
  one-speed-1D
)

set(LIBRARY_NAME cpplib)

foreach(example ${EXAMPLES})
  set(example_name ${example})
  add_executable(${example_name} ${example}.cpp)
  target_link_libraries(${example_name} cpplib)
endforeach()
