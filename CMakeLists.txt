cmake_minimum_required(VERSION 3.22)
project(ParallalSuffixArray)

#set CXX standard to 17



#set(CMAKE_C_COMPILER "clang")
#set(CMAKE_CXX_COMPILER "clang++")

# add flag -std=c++17 to cxx compiler
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

message(STATUS "CMAKE_C_COMPILER: ${CMAKE_C_COMPILER}")
message(STATUS "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")


#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp=libomp")
#find_package(OpenMP)
#if (OPENMP_FOUND)
#    message(STATUS "OpenMP found")
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

#include(FetchContent)
#FetchContent_Declare(
#        parlaylib
#        GIT_REPOSITORY https://github.com/cmuparlay/parlaylib.git
#        GIT_TAG master
#)
#FetchContent_MakeAvailable(parlaylib)
add_executable(build_sa src/main.cpp)

#SET(BUILD_DIVSUFSORT64 ON CACHE BOOL "Build libdivsufsort64")
#add_subdirectory(libdivsufsort)
#add_dependencies(sa divsufsort64)
#target_link_libraries(sa divsufsort64)
#target_include_directories(sa PRIVATE "${CMAKE_BINARY_DIR}/libdivsufsort/include")


set(parlaylib_SOURCE_DIR "/home/alex/Programming/Tools/parlaylib")
target_include_directories(build_sa PRIVATE ${parlaylib_SOURCE_DIR} ${parlaylib_SOURCE_DIR}/include/)

