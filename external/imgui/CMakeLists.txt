cmake_minimum_required(VERSION 2.8.12)
project(imgui)

add_definitions(-DIMGUI_IMPL_OPENGL_LOADER_GLAD)

set(IMGUI_DIR ${CMAKE_CURRENT_SOURCE_DIR}/imgui/)
include_directories(${IMGUI_DIR} ${IMGUI_DIR}/backends ..)

file(GLOB sources ${IMGUI_DIR}/*.cpp ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp ${IMGUI_DIR}/backends/imgui_impl_opengl3.cpp)
add_library(imgui STATIC ${sources})
target_link_libraries(${PROJECT_NAME} glad glfw)
