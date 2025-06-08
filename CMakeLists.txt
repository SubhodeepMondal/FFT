cmake_minimum_required(VERSION 3.20)
project(FFT)


set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

link_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/src/composite_radix_dft.c)