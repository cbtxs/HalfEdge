
add_executable(test_halfedgemesh test_halfedgemesh.cpp)
target_link_libraries(test_halfedgemesh  ${VTK_LIBRARIES})

add_executable(test_pmr test_pmr.cpp)
add_executable(test_chunkarray test_chunkarray.cpp)
add_executable(test_data_container test_data_container.cpp)

add_executable(test_eigen test_eigen.cpp)
target_link_libraries(test_eigen  OpenMP::OpenMP_CXX Eigen3::Eigen)
#target_link_libraries(test_eigen  Eigen3::Eigen)

add_executable(test_plot test_plot.cpp)
target_link_libraries(test_plot ${CAIROMM_LIBRARIES})
target_link_libraries(test_halfedgemesh ${CAIROMM_LIBRARIES})

add_executable(test_simd test_simd.cpp)
target_link_libraries(test_simd ${BLAS_LIBRARIES} OpenMP::OpenMP_CXX Eigen3::Eigen)

add_executable(test_gl test_gl.cpp)
target_link_libraries(test_gl ${GLFW_LIBRARIES} ${OPENGL_LIBRARIES} glad)

add_executable(test_gl_tri test_gl_tri.cpp)
target_link_libraries(test_gl_tri ${GLFW_LIBRARIES} ${OPENGL_LIBRARIES} glad)

add_executable(test_irregulararray2d test_irregulararray2d.cpp)
