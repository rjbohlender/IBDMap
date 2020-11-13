cmake_minimum_required(VERSION 3.0)
project(carvaIBD)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(third_party/SugarPP/)
add_subdirectory(third_party/fmt/)

set(Boost_USE_MULTITHREADED TRUE)

find_package(Armadillo REQUIRED)
find_package(Threads REQUIRED)
find_package(Boost 1.66.0 COMPONENTS iostreams REQUIRED)

if (UNIX)
    find_package(LAPACK REQUIRED)
    find_package(BLAS REQUIRED)
endif()

set(PROJECT_SUPPORT_FILES
    src/jointhreads.hpp
    src/parser.hpp
    src/split.hpp
    src/isgzipped.hpp
    src/permutation.cpp src/permutation.hpp
    src/statistic.cpp src/statistic.hpp
    src/threadpool.hpp
    src/parameters.hpp
    src/threadsafequeue.hpp
    src/reporter.hpp
    src/geneticmap.hpp
    src/indexsort.hpp
    src/breakpoint.hpp
    src/indexer.hpp
    src/glm.hpp
    third_party/stocc/stocc.h
    third_party/stocc/sfmt.cpp third_party/stocc/sfmt.h
    third_party/stocc/stoc3.cpp
    third_party/stocc/fnchyppr.cpp
    third_party/stocc/wnchyppr.cpp
    third_party/stocc/stoc1.cpp
    third_party/stocc/stoc2.cpp
    third_party/stocc/randomc.h
    third_party/stocc/mersenne.cpp
    third_party/stocc/erfres.cpp
    third_party/stocc/userintf.cpp
    link/binomial.cpp
    link/binomial.hpp
    link/family.cpp
    link/family.hpp
    link/gaussian.cpp
    link/gaussian.hpp
    src/info.cpp
    src/info.hpp
    src/math.hpp
    src/math.cpp
    src/parser.cpp
    src/source.cpp
    src/source.hpp)

include_directories("${CMAKE_CURRENT_LIST_DIR}/third_party")
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${LAPACK_INCLUDE_DIRS})
include_directories(${BLAS_INCLUDE_DIRS})

add_executable(carvaIBD main.cpp ${PROJECT_SUPPORT_FILES})

target_link_libraries(carvaIBD fmt::fmt-header-only)
target_link_libraries(carvaIBD pthread)
target_link_libraries(carvaIBD ${LAPACK_LIBRARIES})
target_link_libraries(carvaIBD ${BLAS_LIBRARIES})
target_link_libraries(carvaIBD ${ARMADILLO_LIBRARIES})
target_link_libraries(carvaIBD ${Boost_LIBRARIES})

target_compile_definitions(carvaIBD PRIVATE MAXCOLORS=10000000)
