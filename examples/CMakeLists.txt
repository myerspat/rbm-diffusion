add_executable(examples
  one-speed-1D.cpp
)

target_include_directories(examples PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/
)
