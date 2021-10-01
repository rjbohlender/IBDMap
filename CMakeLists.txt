cmake_minimum_required(VERSION 3.0)
project(carvaIBD)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(third_party/SugarPP/)
add_subdirectory(third_party/fmt/)

set(Boost_USE_MULTITHREADED TRUE)
set(Boost_USE_STATIC_LIBS   ON)

LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")

find_package(Armadillo REQUIRED)
find_package(Threads REQUIRED)
find_package(Boost 1.66.0 COMPONENTS iostreams REQUIRED)
find_package(Gperftools)

SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a") # Find static ZLIB
find_package(ZLIB REQUIRED)
UNSET (CMAKE_FIND_LIBRARY_SUFFIXES)

if (UNIX)
    set(BLA_VENDOR OpenBLAS)
    set(BLA_STATIC ON)
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
        src/source.hpp
        src/inputvalidator.cpp
        src/inputvalidator.hpp)

include_directories("${CMAKE_CURRENT_LIST_DIR}/third_party")
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${LAPACK_INCLUDE_DIRS})
include_directories("${CMAKE_CURRENT_LIST_DIR}/third_party/sse2neon")
include_directories("${CMAKE_CURRENT_LIST_DIR}/third_party/zstr/src")
include_directories(${BLAS_INCLUDE_DIRS})

add_executable(carvaIBD main.cpp ${PROJECT_SUPPORT_FILES})

#target_link_libraries(carvaIBD -static)
target_link_libraries(carvaIBD fmt::fmt-header-only)
target_link_libraries(carvaIBD ${BLAS_LIBRARIES})
target_link_libraries(carvaIBD ${LAPACK_LIBRARIES})
target_link_libraries(carvaIBD ${ARMADILLO_LIBRARIES})
target_link_libraries(carvaIBD ${Boost_LIBRARIES})
target_link_libraries(carvaIBD ${ZLIB_LIBRARIES})
if(Gperftools_FOUND)
    target_link_libraries(carvaIBD ${Gperftools_LIBRARIES})
endif()
target_link_libraries(carvaIBD libpthread.a) # Libpthread is still getting dynamically linked too, but I think that's OK
target_link_libraries(carvaIBD libgfortran.a) # Necessary for static linking BLAS
target_link_libraries(carvaIBD libquadmath.a) # Libgfortran seems to need this

target_compile_definitions(carvaIBD PRIVATE MAXCOLORS=10000000)
