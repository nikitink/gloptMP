cmake_minimum_required(VERSION 3.9)
project(gloptMP)

find_package(inmost REQUIRED)
find_package(OpenMP REQUIRED)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)

link_directories(${INMOST_LIBRARY_DIRS})
include_directories(${INMOST_INCLUDE_DIRS})
add_definitions(${INMOST_DEFINITIONS})

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")

set(SOURCES hilbert.cpp gloptMP.cpp main.cpp)
add_executable(gloptMP ${SOURCES})
target_link_libraries(gloptMP PRIVATE OpenMP::OpenMP_CXX)
