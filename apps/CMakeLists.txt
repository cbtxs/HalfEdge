
add_library(cutmesh MODULE cut_mesh.cpp)
# 静态链接 C++ 标准库
set_target_properties(cutmesh PROPERTIES LINK_FLAGS "-static-libstdc++")

add_executable(cut_mesh cut_mesh.cpp)
