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

list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)


if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/OpenVolumeMesh")
  add_subdirectory(OpenVolumeMesh)
  message("Building with OpenVolumeMesh")
else()
  message("Not Building with OpenVolumeMesh")
endif()

find_package (OpenVolumeMesh REQUIRED)


add_definitions(-DINCLUDE_TEMPLATES)

set (SOURCES "")
set (INCLUDE_DIRS "")
set (LIBRARIES "")
set (LIBRARY_DIRS "")

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
set (HEXEX_COMPILE_FLAGS "-Wall")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set (HEXEX_COMPILE_FLAGS "${HEXEX_COMPILE_FLAGS} -msse -mfpmath=sse" )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set (HEXEX_COMPILE_FLAGS "${HEXEX_COMPILE_FLAGS} -msse -mfpmath=sse -pedantic -Wextra")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    message (WARNING
        "You are using an Intel compiler which might generate x87 FPU code "
        "that breaks the exact predicates. If you know which compiler flags "
        "ensure that the Intel compiler produces SSE code, please patch "
        "the CMakeLists.txt and inform the author <lyon@cs.rwth-aachen.de>.")
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

set_target_properties (HexEx
    PROPERTIES
    COMPILE_FLAGS "${HEXEX_COMPILE_FLAGS}"
    DEFINE_SYMBOLS "-DHEXEX_EXPORT_SYMBOLS"
    CXX_STANDARD 11
    )

target_include_directories(HexEx
        INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        )

set_target_properties (HexExStatic
    PROPERTIES
    COMPILE_FLAGS "${HEXEX_COMPILE_FLAGS}"
    DEFINE_SYMBOLS "-DHEXEX_EXPORT_SYMBOLS"
    CXX_STANDARD 11
    )
set_property(TARGET HexExStatic PROPERTY POSITION_INDEPENDENT_CODE ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

#
# Fake successful finder run if compiling as a dependent project.
#
if(NOT "${PROJECT_NAME}" MATCHES "HexEx")
    set (HEXEX_FOUND true PARENT_SCOPE)
    set (HEXEX_LIBRARIES HexEx PARENT_SCOPE)
    set (HEXEX_LIBRARY HexEx PARENT_SCOPE)
    set (HEXEX_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src/" PARENT_SCOPE)
    set (HEXEX_LIBRARY_DIR "${CMAKE_CURRENT_BINARY_DIR}" CACHE PATH "The directory where the HexEx libraries can be found.")
endif()

if("${PROJECT_NAME}" MATCHES "HexEx")
   add_subdirectory(demo/cmdline_tool)
   add_subdirectory(demo/minimum_example)
endif()

#if ( NOT DEFINED HEXEX_BUILD_UNIT_TESTS)
#  set( HEXEX_BUILD_UNIT_TESTS false CACHE BOOL "Enable or disable unit test builds in HexEx." )
#endif()
#add_subdirectory (unittests)

set(HEXEX_BUILD_UNIT_TESTS false CACHE BOOL "Whether to build the unit tests.")
if (HEXEX_BUILD_UNIT_TESTS)
    enable_testing()
    find_package(GoogleTest REQUIRED)
    add_subdirectory(tests)
endif()


#find_package(Doxygen)
#if (DOXYGEN_FOUND)
#    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)
#    add_custom_target(qex_doc
#        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
#        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
#        COMMENT "Generating Doxygen documentation" VERBATIM)
#endif()

list(REMOVE_AT CMAKE_MODULE_PATH -1)
