cmake_minimum_required(VERSION 3.16)

project(Math3)
set(CMAKE_CXX_STANDARD 17)

# Adiciona os diretórios onde estão os headers
include_directories(
    ${CMAKE_SOURCE_DIR}/includes        # glad/glad.h, glfw, stb_image etc.
    ${CMAKE_SOURCE_DIR}/stb_image
)

link_directories(${PROJECT_SOURCE_DIR}/lib)
include_directories(${PROJECT_SOURCE_DIR}/includes)

# Cria o executável com todos os arquivos fontes necessários
add_executable(Math3
    Math3.cpp
    glad.c
    stb_image/stb_image.cpp
)

target_link_libraries(Math3 glfw3 opengl32)

