cmake_minimum_required(VERSION 3.0.0)

project(HalfEdgeMesh)

set(HalfEdgeMesh_VERSION_MAJOR 0)
set(HalfEdgeMesh_VERSION_MINOR 0)
set(HalfEdgeMesh_VERSION_PATCH 1)

include_directories(${PROJECT_SOURCE_DIR}/include)
set(HalfEdgeMesh_LIBRARIES HEMesh)
message(STATUS "${PROJECT_SOURCE_DIR}")

# requires at least C++14
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wl,--exclude-libs,ALL")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif()
message("CMAKE_BUILD_TYPE : " ${CMAKE_BUILD_TYPE})

add_subdirectory(test)
