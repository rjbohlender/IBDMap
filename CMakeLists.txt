cmake_minimum_required(VERSION 3.14)
project(carvaIBD)

set(CMAKE_CXX_STANDARD 17)

find_package(Armadillo REQUIRED)
find_package(Threads REQUIRED)
find_package(Boost 1.66.0 COMPONENTS program_options iostreams REQUIRED)
find_package(absl REQUIRED)

add_subdirectory(third_party/abseil-cpp)

# if(APPLE)
add_definitions(-DLOWMEM)
# endif(APPLE)

set(PROJECT_SUPPORT_FILES
    src/jointhreads.hpp src/parser.hpp src/split.hpp src/isgzipped.hpp
    src/permutation.cpp src/permutation.hpp third_party/stocc/stocc.h
    third_party/stocc/sfmt.cpp third_party/stocc/sfmt.h third_party/stocc/stoc3.cpp
    third_party/stocc/fnchyppr.cpp src/statistic.cpp src/statistic.hpp src/threadpool.hpp src/parameters.hpp src/threadsafequeue.hpp
    third_party/stocc/wnchyppr.cpp third_party/stocc/stoc1.cpp third_party/stocc/stoc2.cpp third_party/stocc/randomc.h third_party/stocc/mersenne.cpp
    third_party/stocc/erfres.cpp third_party/stocc/userintf.cpp src/reporter.hpp src/geneticmap.hpp)

include_directories("${CMAKE_CURRENT_LIST_DIR}/third_party")
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

add_executable(carvaIBD main.cpp ${PROJECT_SUPPORT_FILES})

target_link_libraries(carvaIBD ${ARMADILLO_LIBRARIES})
target_link_libraries(carvaIBD ${Boost_LIBRARIES})
target_link_libraries(carvaIBD absl::hash absl::synchronization absl::flat_hash_map)
