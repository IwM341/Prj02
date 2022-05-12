cmake_minimum_required(VERSION 3.5)

project(src LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fopenmp")
include_directories( "D:\\Important\\articles\\DMFramework\\utils" )

set(TEST FALSE)
if(${TEST})
    add_executable(src testvel.cpp)
    add_definitions(-DTEST1)
else()
    add_definitions(-DDELTAMASS)
    add_executable(src main.cpp)
endif()
