#
# Copyright 2016 Computer Graphics Group, RWTH Aachen University
# Author: Hans-Christian Ebke <ebke@cs.rwth-aachen.de>
#         Max Lyon <lyon@cs.rwth-aachen.de
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

cmake_minimum_required(VERSION 3.4)

if(NOT TARGET OpenVolumeMesh::OpenVolumeMesh)
    find_package (OpenVolumeMesh REQUIRED)
endif()

add_executable (minimum_example main.cc)
target_link_libraries(minimum_example OpenVolumeMesh::OpenVolumeMesh HexEx)

# disable OpenVolumeMesh debug output:
target_compile_definitions(minimum_example PUBLIC NDEBUG)


target_compile_options(minimum_example PUBLIC
    $<$<CXX_COMPILER_ID:Clang>:-Wno-long-long>
    $<$<CXX_COMPILER_ID:GNU>:-pedantic -Wno-long-long>
    )

