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
    message("Configuring HexEx as standalone project...")
else()
    set(HEXEX_STANDALONE_BUILD FALSE)
    message("Configuring HexEx inside another cmake project...")
endif()

if(HEXEX_STANDALONE_BUILD)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/bin")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/lib")
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/lib")
endif()

#list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)


if(NOT TARGET OpenVolumeMesh::OpenVolumeMesh AND EXISTS "${CMAKE_CURRENT_LIST_DIR}/OpenVolumeMesh")
  add_subdirectory(OpenVolumeMesh)
  message("Building with local OpenVolumeMesh")
else()
  message("Not building with local OpenVolumeMesh")
endif()

if (NOT TARGET OpenVolumeMesh::OpenVolumeMesh)
 find_package (OpenVolumeMesh REQUIRED)
endif()


set(SOURCE_FILES
    src/HexEx/Dart.cc
    src/HexEx/DerivedExactPredicates.cc
    src/HexEx/Direction.cc
    src/HexEx/FileAccessor.cc
    src/HexEx/GridIsomorphism.cc
    src/HexEx/HexEx.cc
    src/HexEx/HexExtractor.cc
    src/HexEx/HPort.cc
    src/HexEx/MergePosition.cc
    src/HexEx/Utils.cc
    src/HexEx/predicates.c)

add_library (HexEx ${SOURCE_FILES})
add_library (HexEx::HexEx ALIAS HexEx)
target_link_libraries (HexEx PUBLIC OpenVolumeMesh::OpenVolumeMesh)

include(GenerateExportHeader)
generate_export_header(HexEx
    BASE_NAME HEXEX
    EXPORT_FILE_NAME "src/HexEx/Config/Export.hh"
    DEPRECATED_MACRO_NAME "CMAKE_HEXEX_DEPRECATED"
    NO_DEPRECATED_MACRO_NAME "CMAKE_HEXEX_NO_DEPRECATED"
)

configure_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/src/HexEx/Config/Version.hh.in"
    "${CMAKE_CURRENT_BINARY_DIR}/src/HexEx/Config/Version.hh"
)


target_include_directories(HexEx
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/HexEx> # get rid of this one some time
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/src>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

set_target_properties (HexEx PROPERTIES
    VERSION ${HexEx_VERSION_MAJOR}.${HexEx_VERSION_MINOR}
    SOVERSION ${HexEx_VERSION_MAJOR}.${HexEx_VERSION_MINOR}
    POSITION_INDEPENDENT_CODE ON
    AUTOMOC OFF
    AUTOUIC OFF
    # TODO: turn this on once we have HEXEX_EXPORT macros in all required places:
    #CXX_VISIBILITY_PRESET hidden
    )


target_compile_definitions(HexEx
    PUBLIC
        "INCLUDE_TEMPLATES"
        "$<$<CXX_COMPILER_ID:MSVC>:NOMINMAX>"
        "$<$<CXX_COMPILER_ID:MSVC>:_USE_MATH_DEFINES>"
    )

set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/src/HexEx/predicates.c
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

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES ".*Clang")
    if (NOT ("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "arm64") )
       target_compile_options(HexEx PUBLIC -msse -mfpmath=sse)
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options(HexEx PUBLIC -msse -mfpmath=sse)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    # On 64-bit platforms, MSVC is supposed to only generate x87 code for
    #   'long double' types, not for doubles.
    # For 32-bit platforms, there is code in exactinit() that should
    #   make sure that no x87 extended internal precision is used.
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


include(CMakePackageConfigHelpers)
include(GNUInstallDirs)


set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/HexEx)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/HexExConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/HexExConfig.cmake"
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
    )

write_basic_package_version_file(
    HexExConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
    )

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/HexExConfigVersion.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/HexExConfig.cmake"
    DESTINATION ${INSTALL_CONFIGDIR})


install(DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/src/HexEx/Config"
    DESTINATION include/HexEx)

# Install Header Files
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/
    DESTINATION include
    FILES_MATCHING PATTERN "*.hh"
    )

install(TARGETS HexEx
    EXPORT HexExTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(EXPORT HexExTargets
    FILE HexExTargets.cmake
    NAMESPACE HexEx::
    DESTINATION ${INSTALL_CONFIGDIR}
)

export(EXPORT HexExTargets
    NAMESPACE HexEx::)



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

#
# Fake successful finder run if compiling as a dependent project.
#
if(NOT "${CMAKE_PROJECT_NAME}" MATCHES "HexEx")
    set (HEXEX_FOUND true PARENT_SCOPE)
    set (HEXEX_LIBRARIES HexEx PARENT_SCOPE)
    set (HEXEX_LIBRARY HexEx PARENT_SCOPE)
    set (HEXEX_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src/" PARENT_SCOPE)
    set (HEXEX_LIBRARY_DIR "${CMAKE_CURRENT_BINARY_DIR}" CACHE PATH "The directory where the HexEx libraries can be found.")
endif()

