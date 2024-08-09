set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

add_library(cutmesh MODULE cut_mesh.cpp)

# 静态链接 C++ 标准库
set_target_properties(cutmesh PROPERTIES LINK_FLAGS "-static-libgcc -static-libstdc++")

# 为模块添加-fPIC
set_target_properties(cutmesh PROPERTIES POSITION_INDEPENDENT_CODE ON)


add_library(cutmesh_new MODULE cut_mesh.cpp)

# 静态链接 C++ 标准库
set_target_properties(cutmesh_new PROPERTIES LINK_FLAGS "-static-libgcc -static-libstdc++")

# 为模块添加-fPIC
set_target_properties(cutmesh_new PROPERTIES POSITION_INDEPENDENT_CODE ON)



add_executable(cut_mesh cut_mesh.cpp )
target_link_libraries(cut_mesh ${CAIROMM_LIBRARIES})

add_executable(cut_mesh_new cut_mesh_new.cpp )
target_link_libraries(cut_mesh_new ${CAIROMM_LIBRARIES})
