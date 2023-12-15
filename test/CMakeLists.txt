
add_executable(test_halfedgemesh test_halfedgemesh.cpp)
add_executable(test_pmr test_pmr.cpp)
add_executable(test_chunkarray test_chunkarray.cpp)
add_executable(test_data_container test_data_container.cpp)

add_executable(test_plot test_plot.cpp)
target_link_libraries(test_plot ${CAIROMM_LIBRARIES})
