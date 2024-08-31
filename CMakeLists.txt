cmake_minimum_required(VERSION 3.27.2)

project(matlibr)

set(CMAKE_C_COMPILER g++)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set(CMAKE_PREFIX_PATH ${CMAKE_SOURCE_DIR}/third_party)
set(PROJECT_ROOT_DIR ${CMAKE_SOURCE_DIR})
set(PROJECT_LIBRARY_NAME matlibr)

add_subdirectory(src)

enable_testing()

add_subdirectory(tests)
