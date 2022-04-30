cmake_minimum_required(VERSION 3.5)

project(src LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fopenmp")
include_directories( "D:\\Important\\articles\\DMFramework\\utils" )
add_executable(src main.cpp)
