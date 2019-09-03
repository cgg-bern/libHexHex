/*
 * Copyright 2016 Computer Graphics Group, RWTH Aachen University
 * Author: Max Lyon <lyon@cs.rwth-aachen.de>
 *
 * This file is part of HexEx.
 *
 * HexEx is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * HexEx is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HexEx.  If not, see <http://www.gnu.org/licenses/>.
 */


#pragma once

#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>
#include <OpenVolumeMesh/Core/BaseEntities.hh>
#include <OpenVolumeMesh/Core/PropertyPtr.hh>
#include <OpenVolumeMesh/Core/PropertyDefines.hh>
#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralGeometryKernel.hh>
#include <OpenVolumeMesh/Core/GeometryKernel.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMeshTopologyKernel.hh>

#include "GridIsomorphism.hh"
#include "Matrix4x4T.hh"

namespace HexEx
{

using Vec2d = OpenVolumeMesh::Geometry::Vec2d;
using Vec3d = OpenVolumeMesh::Geometry::Vec3d;
using Vec4d = OpenVolumeMesh::Geometry::Vec4d;
using Vec3i = OpenVolumeMesh::Geometry::Vec3i;
using Parameter = Vec3d;
using Position = Vec3d;

using TetrahedralMesh = OpenVolumeMesh::TetrahedralGeometryKernel<Vec3d, OpenVolumeMesh::TetrahedralMeshTopologyKernel>;
using PolyhedralMesh = OpenVolumeMesh::GeometryKernel<Vec3d, OpenVolumeMesh::TopologyKernel>;
using HexahedralMesh = OpenVolumeMesh::GeometryKernel<Vec3d, OpenVolumeMesh::HexahedralMeshTopologyKernel>;


using VertexHandle   = OpenVolumeMesh::VertexHandle;
using HalfEdgeHandle = OpenVolumeMesh::HalfEdgeHandle;
using EdgeHandle     = OpenVolumeMesh::EdgeHandle;
using HalfFaceHandle = OpenVolumeMesh::HalfFaceHandle;
using FaceHandle     = OpenVolumeMesh::FaceHandle;
using CellHandle     = OpenVolumeMesh::CellHandle;

using Edge = OpenVolumeMesh::OpenVolumeMeshEdge;
using Face = OpenVolumeMesh::OpenVolumeMeshFace;
using Cell = OpenVolumeMesh::OpenVolumeMeshCell;

template <typename T>
using VertexProperty = OpenVolumeMesh::VertexPropertyT<T>;
template <typename T>
using EdgeProperty = OpenVolumeMesh::EdgePropertyT<T>;
template <typename T>
using HalfEdgeProperty = OpenVolumeMesh::HalfEdgePropertyT<T>;
template <typename T>
using HalfFaceProperty = OpenVolumeMesh::HalfFacePropertyT<T>;
template <typename T>
using FaceProperty = OpenVolumeMesh::FacePropertyT<T>;
template <typename T>
using CellProperty = OpenVolumeMesh::CellPropertyT<T>;

template <typename T>
using VertexMapProp = std::map<VertexHandle, T>;

template <typename T>
using PerCellVertexProperty = OpenVolumeMesh::CellPropertyT<VertexMapProp<T> >;

//using Transition = ACG::Matrix4x4dd;
using Matrix4x4dd = HexEx::Matrix4x4T<double>;

using Transition = GridIsomorphism;


using CellTransitionPair = std::pair<CellHandle, Transition>;

class HPortHandle;
using PortTransitionPair = std::pair<HPortHandle, Transition>;

}
