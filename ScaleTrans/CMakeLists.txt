cmake_minimum_required(VERSION 3.1)
project(ScaleTrans)

add_executable(${PROJECT_NAME}_bin main.cpp DeformationGradient2d.cpp DeformationGradient2d.h DeformableMesh.h TriangleMesh.h TriangleMesh.cpp)
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::opengl igl::opengl_glfw igl::embree)
