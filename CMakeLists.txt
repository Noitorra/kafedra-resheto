cmake_minimum_required(VERSION 2.8)

project(kafedra-resheto)

# Require MPI for this project:
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin/")

ADD_SUBDIRECTORY(src)
