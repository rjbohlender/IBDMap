cmake_minimum_required(VERSION 3.0) # A newer version of cmake will let us statically link MKL
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

set(CMAKE_SHARED_MODULE_PREFIX "")

set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined" )

if (UNIX)
    find_package(LAPACK REQUIRED)
    find_package(BLAS REQUIRED)
    if (NOT APPLE)
        set(BLA_STATIC ON)
    endif()
endif()

set(PROJECT_SUPPORT_FILES
    src/jointhreads.hpp
    src/parser.hpp
    src/split.hpp
    src/isgzipped.hpp
    src/statistic.cpp src/statistic.hpp
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
    src/inputvalidator.hpp src/phenotypes.cpp src/phenotypes.hpp)

set(IBDLIB_SOURCE_FILES
        src/ibdlib.hpp
        src/ibdlib.cpp
        src/split.hpp
        src/geneticmap.hpp
        src/inputvalidator.hpp
        src/inputvalidator.cpp src/ibdlibclasses.hpp)

include_directories("${CMAKE_CURRENT_LIST_DIR}/third_party")
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${Python3_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${LAPACK_INCLUDE_DIRS})
include_directories(${BLAS_INCLUDE_DIRS})
include_directories("${CMAKE_CURRENT_LIST_DIR}/third_party/sse2neon")
include_directories("${CMAKE_CURRENT_LIST_DIR}/third_party/zstr/src")

add_library(ibdlib MODULE ${IBDLIB_SOURCE_FILES})
target_link_libraries(ibdlib ${Python3_LIBRARIES} ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} fmt::fmt-header-only ${ARMADILLO_LIBRARIES})

add_executable(carvaIBD main.cpp ${PROJECT_SUPPORT_FILES})

target_link_libraries(carvaIBD fmt::fmt-header-only)
target_link_libraries(carvaIBD pthread)
target_link_libraries(carvaIBD ${LAPACK_LIBRARIES})
target_link_libraries(carvaIBD ${BLAS_LIBRARIES})
target_link_libraries(carvaIBD ${ARMADILLO_LIBRARIES})
target_link_libraries(carvaIBD ${Python3_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(carvaIBD ${ZLIB_LIBRARIES})
if(Gperftools_FOUND)
  target_link_libraries(carvaIBD ${Gperftools_LIBRARIES})
endif()

target_compile_definitions(carvaIBD PRIVATE MAXCOLORS=10000000)
