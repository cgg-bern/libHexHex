#pragma once

#include <HexHex/Commons/MeshElement.hh>
#include <HexHex/Utils/Typedefs.hh>

namespace HexHex
{

// Forward declaration
class Direction;


HalfFaceHandle findOppositeNextHalfFace (
    const TetMeshCache& mesh,
    const CellHandle ch, const std::array<VertexHandle,4>& cell_vhs, const std::array<Parameter,4>& cell_params,
    const HalfFaceHandle hfh, const Parameter& u,
    const Direction d1, const Direction d2, MeshElement& traceThrough, bool computeTraceThrough = false);

HalfFaceHandle findBladeNextHalfFace (
    const TetMeshCache& mesh,
    const CellHandle ch, const std::array<VertexHandle,4>& cell_vhs, const std::array<Parameter,4>& cell_params,
    const MeshElement& generator, const MeshElement& holder,
    const HalfFaceHandle hfh, const Parameter& u,
    const Direction d1, const Direction d2);

}

