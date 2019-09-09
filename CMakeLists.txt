cmake_minimum_required(VERSION 3.14)
project(carvaIBD)

set(CMAKE_CXX_STANDARD 17)

find_package(Armadillo REQUIRED)
find_package(Threads REQUIRED)
find_package(Boost 1.66.0 COMPONENTS program_options iostreams REQUIRED)

set(PROJECT_SUPPORT_FILES
    src/jointhreads.hpp src/parser.cpp src/parser.hpp src/split.hpp)

include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})


add_executable(carvaIBD main.cpp ${PROJECT_SUPPORT_FILES})