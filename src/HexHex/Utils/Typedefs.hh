#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Core/BaseEntities.hh>
#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralGeometryKernel.hh>
#include <OpenVolumeMesh/Core/GeometryKernel.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMeshTopologyKernel.hh>

//#include <absl/hash/hash.h>

namespace HexHex
{
class HexExtractor;
class Direction;
class CellIGM;
class Transition;
class MeshElement;
class Parametrization;
struct GlobalPropellerHandle;
struct TetMeshCache;
struct PLMesh;

using int32 = int32_t;
using uint32 = uint32_t;
using uint = unsigned int;

namespace OVM = OpenVolumeMesh;

using Vec2d = OVM::Vec2d;
using Vec3d = OVM::Vec3d;
using Vec3i = OVM::Vec3i;
using Parameter = Vec3d;
using Position = Vec3d;

template <typename VecT>
using TetrahedralMeshT = OpenVolumeMesh::TetrahedralGeometryKernel<VecT, OpenVolumeMesh::TetrahedralMeshTopologyKernel>;
template <typename VecT>
using PolyhedralMeshT = OpenVolumeMesh::GeometryKernel<VecT, OpenVolumeMesh::TopologyKernel>;
template<typename VecT>
using HexahedralMeshT = OpenVolumeMesh::GeometryKernel<VecT, OpenVolumeMesh::HexahedralMeshTopologyKernel>;

using TetrahedralMesh = TetrahedralMeshT<Vec3d>;
using PolyhedralMesh = PolyhedralMeshT<Vec3d>;
using HexahedralMesh = HexahedralMeshT<Vec3d>;

using VertexHandle   = OpenVolumeMesh::VertexHandle;
using HalfEdgeHandle = OpenVolumeMesh::HalfEdgeHandle;
using EdgeHandle     = OpenVolumeMesh::EdgeHandle;
using HalfFaceHandle = OpenVolumeMesh::HalfFaceHandle;
using FaceHandle     = OpenVolumeMesh::FaceHandle;
using CellHandle     = OpenVolumeMesh::CellHandle;

using HFParam = OVM::HalfFacePropertyT<HexHex::Vec3d>;

} // namespace HexHex

