cmake_minimum_required(VERSION 3.11) # A newer version of cmake will let us statically link MKL
set(CMAKE_OSX_DEPLOYMENT_TARGET 12.3) # Otherwise it will pull an older SDK, which breaks linking against homebrew packages
project(carvaIBD)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(third_party/SugarPP/)
add_subdirectory(third_party/fmt/ EXCLUDE_FROM_ALL)

set(Boost_USE_MULTITHREADED TRUE)

LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(Armadillo REQUIRED)
find_package(Threads REQUIRED)
find_package(Boost 1.66.0 COMPONENTS iostreams python numpy REQUIRED)
find_package(Gperftools)
find_package(ZLIB REQUIRED)

find_package(benchmark)

set(CMAKE_SHARED_MODULE_PREFIX "")

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined" )
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -DNDEBUG -fno-omit-frame-pointer")

if (UNIX)
    find_package(LAPACK REQUIRED)
    find_package(BLAS REQUIRED)
endif()

set(PROJECT_SUPPORT_FILES
    third_party/stocc/stocc.h
    third_party/stocc/stoc1.cpp
    third_party/stocc/stoc2.cpp
    third_party/stocc/stoc3.cpp
    third_party/stocc/fnchyppr.cpp
    third_party/stocc/wnchyppr.cpp
    third_party/stocc/erfres.cpp
    third_party/stocc/randomc.h
    third_party/stocc/mersenne.cpp
    third_party/stocc/userintf.cpp
    third_party/stocc/sfmt.h
    third_party/stocc/sfmt.cpp
    src/jointhreads.hpp
    src/parser.hpp
    src/split.hpp
    src/iscompressed.hpp
    src/statistic.cpp
    src/statistic.hpp
    src/threadpool.hpp
    src/parameters.hpp
    src/threadsafequeue.hpp
    src/reporter.hpp
    src/geneticmap.hpp
    src/indexsort.hpp
    src/breakpoint.hpp
    src/indexer.hpp
    src/info.cpp
    src/info.hpp
    src/math.hpp
    src/math.cpp
    src/parser.cpp
    src/source.cpp
    src/source.hpp
    src/inputvalidator.cpp
    src/inputvalidator.hpp
    src/phenotypes.cpp
    src/phenotypes.hpp
    src/glm.hpp
    src/permutation.hpp
    src/permutation.cpp
    link/family.hpp
    link/family.cpp
    link/binomial.hpp
    link/binomial.cpp
    link/gaussian.hpp
    link/gaussian.cpp)

set(IBDLIB_SOURCE_FILES
    src/ibdlib.hpp
    src/ibdlib.cpp
    src/split.hpp
    src/geneticmap.hpp
    src/inputvalidator.hpp
    src/inputvalidator.cpp src/ibdlibclasses.hpp
    src/glm.hpp
    src/permutation.hpp
    src/permutation.cpp)

include_directories("${CMAKE_CURRENT_LIST_DIR}/third_party")
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${Python3_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${LAPACK_INCLUDE_DIRS})
include_directories(${BLAS_INCLUDE_DIRS})
include_directories("${CMAKE_CURRENT_LIST_DIR}/third_party/sse2neon")
include_directories("${CMAKE_CURRENT_LIST_DIR}/third_party/zstr/src")
include_directories("${CMAKE_CURRENT_LIST_DIR}/third_party/pcg-cpp-0.98/include")

add_library(ibdlib MODULE ${IBDLIB_SOURCE_FILES})
target_link_libraries(ibdlib ${Python3_LIBRARIES} ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} fmt::fmt-header-only ${ARMADILLO_LIBRARIES})

add_executable(IBDMap main.cpp ${PROJECT_SUPPORT_FILES})

target_compile_definitions(IBDMap PRIVATE ARMA_DONT_USE_WRAPPER=1)

target_link_libraries(IBDMap fmt::fmt-header-only)
target_link_libraries(IBDMap pthread)
target_link_libraries(IBDMap ${LAPACK_LIBRARIES})
target_link_libraries(IBDMap ${BLAS_LIBRARIES})
target_link_libraries(IBDMap ${Python3_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(IBDMap ${ZLIB_LIBRARIES})
if(Gperftools_FOUND)
  target_link_libraries(IBDMap ${Gperftools_LIBRARIES})
endif()

target_compile_definitions(IBDMap PRIVATE MAXCOLORS=10000000)

if(benchmark_FOUND)
    add_executable(ibdBench bench_main.cpp ${PROJECT_SUPPORT_FILES})
    target_compile_definitions(ibdBench PRIVATE ARMA_DONT_USE_WRAPPER=1)

    target_link_libraries(ibdBench benchmark)
    target_link_libraries(ibdBench fmt::fmt-header-only)
    target_link_libraries(ibdBench pthread)
    target_link_libraries(ibdBench ${BLAS_LIBRARIES})
    target_link_libraries(ibdBench ${LAPACK_LIBRARIES})
    target_link_libraries(ibdBench ${Python3_LIBRARIES} ${Boost_LIBRARIES})
    target_link_libraries(ibdBench ${ZLIB_LIBRARIES})
endif()
