# Copyright 2025 Computer Graphics Group, University of Bern
# Copyright 2013 Computer Graphics Group, RWTH Aachen University
# Authors: Hans-Christian Ebke <ebke@cs.rwth-aachen.de>
#          Max Lyon <lyon@cs.rwth-aachen.de>
#          Tobias Kohler <tobias.kohler@unibe.ch>
#          Martin Heistermann <martin.heistermann@unibe.ch>
#
# This file is part of HexHex.
#
# HexHex is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# HexHex is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU General Public License
# along with HexHex.  If not, see <http://www.gnu.org/licenses/>.
#

cmake_minimum_required(VERSION 3.25 FATAL_ERROR)
cmake_policy(VERSION 3.25.0)

project (HexHexSuperproject
    LANGUAGES CXX C
    HOMEPAGE_URL http://github.com/cgg-bern/libHexHex
    )
set(CMAKE_CXX_STANDARD 20)

if( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
    set(HEXHEX_STANDALONE_BUILD TRUE)
    message("Configuring HexHex as standalone project...")
else()
    set(HEXHEX_STANDALONE_BUILD FALSE)
    message("Configuring HexHex inside another cmake project...")
endif()

if(HEXHEX_STANDALONE_BUILD)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/bin")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/lib")
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/lib")
endif()

#set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)

option(HEXHEX_WITH_JEMALLOC "Build with jemalloc, enables precise memory measurements" ${HEXHEX_STANDALONE_BUILD})


include(FetchContent)
set(FETCHCONTENT_UPDATES_DISCONNECTED TRUE)
set(EXTERNAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external")


if (NOT TARGET OpenVolumeMesh::OpenVolumeMesh)
    FetchContent_Declare(ovm
        GIT_REPOSITORY https://www.graphics.rwth-aachen.de:9000/OpenVolumeMesh/OpenVolumeMesh
        GIT_TAG v3.4.1
        SOURCE_DIR "${EXTERNAL_DIR}/OpenVolumeMesh"
        )
    FetchContent_MakeAvailable(ovm)
    message("Downloaded OpenVolumeMesh to ${ovm_SOURCE_DIR}")
endif()

if (NOT TARGET Timekeeper::libTimekeeper)
    FetchContent_Declare(timekeeper
        GIT_REPOSITORY https://github.com/cgg-bern/libTimekeeper
        GIT_TAG main
        SOURCE_DIR "${EXTERNAL_DIR}/timekeeper"
        )
    FetchContent_MakeAvailable(timekeeper)
    message("Downloaded libTimekeeper to ${timekeeper_SOURCE_DIR}")
    #add_subdirectory("${timekeeper_SOURCE_DIR}/src")
endif()

if (NOT TARGET nlohmann::json)
    FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz)
    FetchContent_MakeAvailable(json)
    message("Downloaded nlohmann::json to ${json_SOURCE_DIR}")
endif()

if (NOT TARGET absl::hash)
    FetchContent_Declare(
        abseil
        EXCLUDE_FROM_ALL
        GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git
        SOURCE_DIR "${EXTERNAL_DIR}/abseil"
    )
    FetchContent_MakeAvailable(abseil)
    message("Downloaded abseil to ${abseil_SOURCE_DIR}")
endif()

option(HEXHEX_ENABLE_OPENMP "HexHex: try to build with OpenMP" ON)
if (HEXHEX_ENABLE_OPENMP)
    if(APPLE)
        set(OpenMP_ROOT "/opt/homebrew/opt/libomp")
    endif()
    set(OpenMP_RUNTIME_MSVC "llvm") # required for `#prgama omp atomic write`
    find_package(OpenMP COMPONENTS CXX)
endif()

if(HEXHEX_WITH_JEMALLOC)
    find_package(PkgConfig)
    if (PKG_CONFIG_FOUND)
        pkg_check_modules (JEMALLOC jemalloc)
        pkg_search_module(jemalloc IMPORTED_TARGET jemalloc)
    else()
        message(STATUS "pkgconfig not found, building without jemalloc")
    endif()
endif()

add_subdirectory(src)

set(HEXHEX_BUILD_TOOLS ${PROJECT_IS_TOP_LEVEL} CACHE BOOL "Build HexHex commandline tools")
if(HEXHEX_BUILD_TOOLS)

    if(NOT TARGET CLI11::CLI11)
        FetchContent_Declare(cli11
            GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
            GIT_TAG        v2.5.0
            SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/CLI11"
            )
        FetchContent_MakeAvailable(cli11)
        message("Downloaded CLI11 to ${cli11_SOURCE_DIR}")
    endif()

    add_subdirectory(tools/HexHex)
    add_subdirectory(tools/minimum_example)
    add_subdirectory(tools/experiments)
endif()


#=========================
# Testing
#=========================
set(HEXHEX_BUILD_UNIT_TESTS true CACHE BOOL "Whether to build the unit tests.")
if (HEXHEX_STANDALONE_BUILD AND HEXHEX_BUILD_UNIT_TESTS)
    enable_testing()
    if (BUILD_SHARED_LIBS)
        message(WARNING "HexHex unittests can only be linked against a static libHexHex, however BUILD_SHARED_LIBS is enabled. Skipping tests.")
    else()
        add_subdirectory(tests)
    endif()
endif()


#===================================================================
# Fake successful finder run if compiling as a dependent project.
#===================================================================
if(NOT "${CMAKE_PROJECT_NAME}" MATCHES "HexHex")
    set (HEXHEX_FOUND true PARENT_SCOPE)
    set (HEXHEX_LIBRARIES HexHex::HexHex PARENT_SCOPE)
    set (HEXHEX_LIBRARY HexHex::HexHex PARENT_SCOPE)
endif()
