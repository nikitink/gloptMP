cmake_minimum_required(VERSION 3.12)
project(gloptMP)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")

link_directories(${INMOST_LIBRARY_DIRS})
include_directories(${INMOST_INCLUDE_DIRS})
add_definitions(${INMOST_DEFINITIONS})

include_directories(${PYTHON_INCLUDE_DIRS})

add_library(lib3ph SHARED IMPORTED)
set_property(TARGET lib3ph PROPERTY IMPORTED_LOCATION "/home/nikitink/Work/filter-project/3ph-build/libMULTIPHASELIB.so")
list(APPEND LIB_NAMES lib3ph)

find_package(inmost REQUIRED)
find_package(OpenMP REQUIRED)

find_package(GMSH REQUIRED)
list(APPEND LIB_NAMES GMSH::GMSH)

find_package(avsim REQUIRED)
list(APPEND LIB_NAMES VirtualSuturer)

set(SOURCES hilbert.cpp gmshfrac.cpp gloptMP.cpp main.cpp)
add_executable(gloptMP ${SOURCES})
target_link_libraries(gloptMP PUBLIC OpenMP::OpenMP_CXX ${INMOST_LIBRARIES} ${LIB_NAMES})
