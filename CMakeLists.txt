# Require CMAKE 3.1 of higher
cmake_minimum_required(VERSION 3.10.0)

# Project name
project(rbm-diffusion VERSION 0.1.0)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(LIBRARY_NAME cpplib)

include(CTest)
include_directories(src)

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(examples)

enable_testing()
