cmake_minimum_required(VERSION 3.1)
project(triangulation CXX)

add_subdirectory(external)
include_directories(external)

add_executable(triangulation main.cpp)

target_link_libraries(triangulation PUBLIC  cgalsource  malloccountfiles gmp  ${CMAKE_DL_LIBS})
