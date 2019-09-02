# Copyright 2013 Computer Graphics Group, RWTH Aachen University
# Authors: Hans-Christian Ebke <ebke@cs.rwth-aachen.de>
#          Max Lyon <lyon@cs.rwth-aachen.de>
#
# This file is part of HexEx.
#
# HexEx is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# HexEx is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU General Public License
# along with HexEx.  If not, see <http://www.gnu.org/licenses/>.
#

cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
cmake_policy(VERSION 3.12.0)

project (HexEx
    VERSION 1.0.0
    LANGUAGES CXX C
    HOMEPAGE_URL https://www.graphics.rwth-aachen.de/software/libHexEx
    )


if( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
    set(HEXEX_STANDALONE_BUILD TRUE)
else()
    set(HEXEX_STANDALONE_BUILD FALSE)
endif()

#list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)


if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/OpenVolumeMesh")
  add_subdirectory(OpenVolumeMesh)
  message("Building with local OpenVolumeMesh")
else()
  message("Not building with local OpenVolumeMesh")
endif()

find_package (OpenVolumeMesh REQUIRED)


set(SOURCE_FILES
    src/Dart.cc
    src/DerivedExactPredicates.cc
    src/Direction.cc
    src/FileAccessor.cc
    src/GridIsomorphism.cc
    src/HexEx.cc
    src/HexExtractor.cc
    src/HPort.cc
    src/MergePosition.cc
    src/Utils.cc
    src/predicates.c)

add_library (HexEx ${SOURCE_FILES})
target_link_libraries (HexEx OpenVolumeMesh::OpenVolumeMesh)

target_compile_definitions(HexEx
    PUBLIC
        "INCLUDE_TEMPLATES"
        "$<$<CXX_COMPILER_ID:MSVC>:NOMINMAX>"
        "$<$<CXX_COMPILER_ID:MSVC>:_USE_MATH_DEFINES>"
    )

set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/src/predicates.c
    PROPERTIES
    # make "if MSVC" work
    COMPILE_DEFINITIONS "$<$<CXX_COMPILER_ID:MSVC>:MSVC=1>"
# Shewchuk's predicates.c use some old-style C, we don't need warnings about it:
    COMPILE_FLAGS "$<$<CXX_COMPILER_ID:MSVC>:/wd4131>"
    )

#
# In order for the exact predicates to work the compiler
# must not generate x87 FPU code as this leads to the use
# of extended precision registers which lead to
# wrong results.
#
# As SSE does not have extended precision registers,
# forcing the generation of SSE code ensures that the
# exact predicates produce correct results.
#

target_compile_options(HexEx
    PUBLIC
    $<$<CXX_COMPILER_ID:Clang>:-msse -mfpmath=sse>
    $<$<CXX_COMPILER_ID:AppleClang>:-msse -mfpmath=sse>
    $<$<CXX_COMPILER_ID:GNU>:-msse -mfpmath=sse>
)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    # There is code in exactinit() that should make sure that nor
    # x87 extended internal precision is used.
else ()
    message (WARNING
        "You are using an unknown compiler which might generate x87 FPU code "
        "that breaks the exact predicates. If you know how to detect this compiler "
        "and which flags "
        "ensure that this compiler produces SSE code, please patch "
        "the CMakeLists.txt and inform the author <lyon@cs.rwth-aachen.de>.")
endif ()

set(STL_RANGE_CHECKS false CACHE BOOL "Include STL range checks in debug mode (This option is only used in debug mode.)")
# Add a flag to check stl vectors in debugging mode
if (STL_RANGE_CHECKS)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC" )
endif()


target_compile_options(HexEx
    PUBLIC
    "$<$<CXX_COMPILER_ID:MSVC>:/bigobj>"
    "$<$<CXX_COMPILER_ID:MSVC>:/Zc:__cplusplus>" # Enable updated __cplusplus macro

    PRIVATE
    "$<$<CXX_COMPILER_ID:GNU>:-Wall>"
    "$<$<CXX_COMPILER_ID:Clang>:-Wall;-Wno-c++98-compat;-Wno-c++98-compat-pedantic;-Wno-weak-vtables;-Wno-padded;-Wno-sign-conversion;-Wno-missing-prototypes>"
    )

target_include_directories(HexEx
        INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        )

set_target_properties (HexEx
    PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    POSITION_INDEPENDENT_CODE ON
    # TODO: turn this on once we have HEXEX_EXPORT macros in all required places:
    #CXX_VISIBILITY_PRESET hidden
    )


if(HEXEX_STANDALONE_BUILD)
   add_subdirectory(demo/cmdline_tool)
   add_subdirectory(demo/minimum_example)
endif()


set(HEXEX_BUILD_UNIT_TESTS false CACHE BOOL "Whether to build the unit tests.")
if (HEXEX_BUILD_UNIT_TESTS)
    enable_testing()
    find_package(GTest REQUIRED)
    add_subdirectory(tests)
endif()



#list(REMOVE_AT CMAKE_MODULE_PATH -1)
