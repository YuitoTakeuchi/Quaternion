cmake_minimum_required(VERSION 4.1.2)

set(LIBRARY_NAME "quaternion")

# set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

add_library(${LIBRARY_NAME} INTERFACE include/Quaternion.hpp)

target_compile_features(${LIBRARY_NAME} PRIVATE cxx_std_17)