cmake_minimum_required(VERSION 3.16)

project(namfinder VERSION 0.1.1)
include(FetchContent)


option(ENABLE_AVX "Enable AVX2 support" OFF)

find_package(ZLIB)
find_package(Threads)
find_package(OpenMP)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: RelWithDebInfo Debug Release" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "RelWithDebInfo" "Debug" "Release")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
add_compile_options(-Wall -Wextra -Werror=maybe-uninitialized)

# Obtain version from Git or fall back to PROJECT_VERSION if not building
# from a Git repository
add_custom_target(version
        ${CMAKE_COMMAND}
        -D CONFIGIN="${PROJECT_SOURCE_DIR}/src/version.hpp.in"
        -D CONFIGOUT="${PROJECT_BINARY_DIR}/version.hpp"
        -D DEFAULT_VERSION="${PROJECT_VERSION}"
        -P ${CMAKE_SOURCE_DIR}/GitVersion.cmake
        )

configure_file(
        "${PROJECT_SOURCE_DIR}/src/buildconfig.hpp.in"
        "${PROJECT_BINARY_DIR}/buildconfig.hpp"
)

add_library(salib STATIC ${SOURCES}
  src/refs.cpp
  src/fastq.cpp
  src/cmdline.cpp
  src/index.cpp
  src/indexparameters.cpp
  src/output.cpp
  src/pc.cpp
  src/aln.cpp
  #src/cigar.cpp
  #src/aligner.cpp
  src/nam.cpp
  src/randstrobes.cpp
  src/version.cpp
  src/io.cpp
  ext/xxhash.c
)
target_include_directories(salib PUBLIC src/ ext/ ${PROJECT_BINARY_DIR})
target_link_libraries(salib PUBLIC ZLIB::ZLIB Threads::Threads)
IF(ENABLE_AVX)
  target_compile_options(salib PUBLIC "-mavx2")
ENDIF()
add_dependencies(salib version)

add_executable(namfinder src/main.cpp)
target_link_libraries(namfinder PUBLIC salib)
install(TARGETS namfinder DESTINATION bin)




