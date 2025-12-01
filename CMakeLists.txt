cmake_minimum_required(VERSION 3.10)
project(ConicIntersection)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB SRC_FILES src/*.cpp)

add_executable(ConicIntersection
    main.cpp
    ${SRC_FILES}
)

target_include_directories(ConicIntersection PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen
)