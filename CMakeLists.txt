project(simple_example)
cmake_minimum_required(VERSION 3.14)
set(SOURCE_FILES knap_MPI.c)
set(CMAKE_C_COMPILER mpicc)
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
add_executable(simple_example ${SOURCE_FILES})
target_link_libraries(simple_example ${MPI_C_LIBRARIES})

