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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif()
message("CMAKE_BUILD_TYPE : " ${CMAKE_BUILD_TYPE})

# Eigen
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
if(NOT ${Eigen3_FOUND})
    message("Eigen3 not found")
else()
    message("Found Eigen3!")
endif()

# 查找 Cairo 库
find_package(PkgConfig REQUIRED)
pkg_check_modules(CAIROMM REQUIRED cairomm-1.0)
if (CAIROMM_FOUND)
    message(STATUS "Cairo found. Version: ${CAIROMM_VERSION}")
    include_directories(${CAIROMM_INCLUDE_DIRS})
    link_directories(${CAIROMM_LIBRARY_DIRS})
else()
    message(FATAL_ERROR "Cairo not found. Please install Cairo and try again.")
endif()

# 查找 GTK3 库
pkg_check_modules(GTK3 REQUIRED gtk+-3.0)
if (GTK3_FOUND)
    message(STATUS "GTK3 found. Version: ${GTK3_VERSION}")
    include_directories(${GTK3_INCLUDE_DIRS})
    link_directories(${GTK3_LIBRARY_DIRS})
else()
    message(FATAL_ERROR "GTK not found. Please install GTK and try again.")
endif()

# GLAD 库
include_directories(thirdparty/include)
add_library(glad thirdparty/src/glad.c)

# GLFW 库
pkg_check_modules(GLFW REQUIRED glfw3)
if (GLFW_FOUND)
    message(STATUS "GLFW found. Version: ${GLFW_VERSION}")
    include_directories(${GLFW_INCLUDE_DIRS})
    link_directories(${GLFW_LIBRARY_DIRS})
else()
    message(FATAL_ERROR "GLFW not found. Please install GLFW and try again.")
endif()

# BLAS 库
find_package(BLAS REQUIRED)
if (BLAS_FOUND)
    message(STATUS "BLAS found. Version: ${BLAS_VERSION}")
    include_directories(${BLAS_INCLUDE_DIRS})
    link_directories(${BLAS_LIBRARY_DIRS})
else()
    message(FATAL_ERROR "BLAS not found. Please install BLAS and try again.")
endif()

# OpenMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    message(STATUS "OPENMP found. Version: ${OPENMP_VERSION}")
    include_directories(${OPENMP_INCLUDE_DIRS})
    link_directories(${OPENMP_LIBRARY_DIRS})
else()
    message(FATAL_ERROR "OPENMP not found. Please install OPENMP and try again.")
endif()


# 查找 OpenGL 库
set(OpenGL_GL_PREFERENCE "GLVND")
find_package(OpenGL REQUIRED)
if (OpenGL_FOUND)
    message(STATUS "OpenGL found. Version: ${OPENGL_VERSION}")
    include_directories(${OPENGL_INCLUDE_DIRS})
    link_directories(${OPENGL_LIBRARY_DIRS})
else()
    message(FATAL_ERROR "OpenGL not found. Please install OpenGL and try again.")
endif()

add_subdirectory(test)
add_subdirectory(apps)
