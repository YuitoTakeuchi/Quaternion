cmake_minimum_required(VERSION 4.1.2)

set(LIBRARY_NAME "quaternion")

add_library(${LIBRARY_NAME} INTERFACE)

target_include_directories(${LIBRARY_NAME} INTERFACE include)

target_compile_features(${LIBRARY_NAME} INTERFACE cxx_std_17)