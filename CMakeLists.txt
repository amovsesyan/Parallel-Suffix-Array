project(ParallalSuffixArray)

cmake_minimum_required(VERSION 3.25)


set(CMAKE_C_COMPILER "/opt/homebrew/Cellar/llvm/17.0.6_1/bin/clang")
set(CMAKE_CXX_COMPILER "/opt/homebrew/Cellar/llvm/17.0.6_1/bin/clang++")

message(STATUS "CMAKE_C_COMPILER: ${CMAKE_C_COMPILER}")
message(STATUS "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")


set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp=libomp")
find_package(OpenMP)
if (OPENMP_FOUND)
    message(STATUS "OpenMP found")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
add_executable(sa src/main.cpp)
