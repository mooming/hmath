# CMakeList.txt : CMake project for hmath, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project ("hmath")

set(CMAKE_CXX_FLAGS "-std=c++17")

# Add source to this project's executable.
file(GLOB SRC_FILES "*.h" "*.cpp")
add_executable (hmath ${SRC_FILES})

# TODO: Add tests and install targets if needed.
