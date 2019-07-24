cmake_minimum_required(VERSION 3.14)
project(test_fixed_point)

include_directories(
    include
)

add_executable(
    ${PROJECT_NAME}
    src/main.cpp
)